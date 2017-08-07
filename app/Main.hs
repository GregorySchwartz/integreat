{- integrate
Gregory W. Schwartz

Integrate data from multiple sources to find consistent (or inconsistent)
entities.
-}

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE QuasiQuotes #-}

module Main where

-- Standard
import Control.Monad
import Data.Bool
import qualified Data.Map.Strict as Map
import Data.Maybe
import Data.Monoid
import qualified Data.Set as Set
import System.IO

-- Cabal
import qualified Control.Lens as L
import Control.Monad.Trans
import qualified Data.ByteString.Lazy.Char8 as CL
import qualified Data.Csv as CSV
import qualified Data.Text as T
import qualified Data.Text.IO as T
import qualified Data.Vector as V
import Options.Generic

-- import qualified Foreign.R as R
-- import Language.R.Instance as R
-- import Language.R.QQ

-- Local
import Types
import Utility
import Load
import Edge.Correlation
-- import Edge.Aracne
-- import Alignment.CSRW
import Integrate
import Print

-- | Command line arguments
data Options = Options { dataInput          :: Maybe String
                                           <?> "([STDIN] | FILE) The input file containing the data intensities. Follows the format: dataLevel,dataReplicate,vertex,intensity. dataLevel is the level (the base level for the experiment, like \"proteomic_Myla\" or \"RNA_MyLa\" for instance), dataReplicate is the replicate in that experiment that the entity is from (the name of that data set with the replicate name, like \"RNA_MyLa_1\"), and vertex is the name of the entity (must match those in the vertex-input), and the intensity is the value of this entity in this data set."
                       , vertexInput        :: Maybe String
                                           <?> "([Nothing] | FILE) The input file containing similarities between entities. Follows the format: vertexLevel1,vertexLevel2, vertex1,vertex2,similarity. vertexLevel1 is the level (the base title for the experiment, \"data set\") that vertex1 is from, vertexLevel2 is the level that vertex2 is from, and the similarity is a number representing the similarity between those two entities. If not specified, then the same entity (determined by vertex in data-input) will have a similarity of 1, different entities will have a similarity of 0."
                       , entityDiff         :: Maybe T.Text
                                           <?> "([Nothing] | STRING) When comparing entities that are the same, ignore the text after this separator. Used for comparing phosphorylated positions with another level. For example, if we have a strings ARG29 and ARG29_7 that we want to compare, we want to say that their value is the highest in correlation, so this string would be \"_\""
                       , alignmentMethod    :: Maybe String
                                           <?> "([CosineSimilarity] | RandomWalker | RandomWalkerSim) The method to get integrated vertex similarity between levels. CosineSimilarity uses the cosine similarity of each  vertex in each network compared to the other vertices in  other networks. RandomWalker uses a random walker with restart based network algnment algorithm in order to get similarity. RandomWalkerSim uses a random walker with restart and actually simulates the walker to get a stochastic result."
                       , edgeMethod         :: Maybe String
                                           <?> "([SpearmanCorrelation] | PearsonCorrelation ) The method to use for the edges between entities in the coexpression matrix."
                       , walkerRestart      :: Maybe Double
                                           <?> "([0.25] | PROBABILITY) For the random walker algorithm, the probability of making  a jump to a random vertex. Recommended to be the ratio of  the total number of vertices in the top 99% smallest  subnetworks to the total number of nodes in the reduced  product graph (Jeong, 2015)."
                       , steps              :: Maybe Int
                                           <?> "([100] | STEPS) For the random walker algorithm, the number of steps to take  before stopping."
                       , premade            :: Bool
                                           <?> "([False] | BOOL) Whether the input data (dataInput) is a pre-made network of the format \"[([\"VERTEX\"], [(\"SOURCE\", \"DESTINATION\", WEIGHT)])]\", where VERTEX, SOURCE, and DESTINATION are of type INT starting at 0, in order, and WEIGHT is a DOUBLE representing the weight of the edge between SOURCE and DESTINATION."
                       , test               :: Bool
                                           <?> "([False] | BOOL) Whether the input data from premade is from a test run. If supplied, the output is changed to an accuracy measure. In this case, we get the total rank below the number of permuted vertices divided by the theoretical maximum (so if there were five changed vertices out off 10 and two were rank 8 and 10 while the others were in the top five, we would have (1 - ((3 + 5) / (10 + 9 + 8 + 7 + 6))) as the accuracy."
                       , entityFilter       :: Maybe Int
                                           <?> "([Nothing] | INT) The minimum number of samples an entity must appear in, otherwise the entity is removed from the analysis."
                       , entityFilterStdDev :: Maybe Double
                                           <?> "([Nothing] | DOUBLE) Remove entities that have less than this value for their standard deviation among all samples."
                       , permutations       :: Maybe Int
                                           <?> "([1000] | INT) The number of permutations for cosine similarity permutation test or bootstrap. Right now just does bootstrap and only shows the first comparison if there are multiple comparisons."
                       }
               deriving (Generic)

instance ParseRecord Options

-- | Get all of the required information for integration.
getIntegrationInput
    :: Options
    -> IO (Maybe (Set.Set ID), Maybe UnifiedData, IDMap, IDVec, VertexSimMap, EdgeSimMap, GrMap)
getIntegrationInput opts = do
    let processCsv = snd . either error id

    hPutStrLn stderr "Getting data input."
    dataEntries   <- fmap (processCsv . CSV.decodeByName)
                   . maybe CL.getContents CL.readFile
                   . unHelpful
                   . dataInput
                   $ opts

    let numSamples     = fmap NumSamples . unHelpful . entityFilter $ opts
        stdDevThresh   =
            fmap StdDevThreshold . unHelpful . entityFilterStdDev $ opts
        levels         = (\x -> maybe x (flip filterEntitiesStdDev x) stdDevThresh)
                       . entitiesToLevels
                       . (\x -> maybe x (flip filterEntities x) numSamples)
                       . datasToEntities
                       . V.toList
                       $ dataEntries
        unifiedData    = unifyAllLevels . fmap snd $ levels
        levelNames     = Set.toList . Set.fromList . fmap fst $ levels
        idMap          = getIDMap unifiedData
        idVec          = getIDVec unifiedData
        size           = Size . Map.size . unIDMap $ idMap
        eDiff          = fmap EntityDiff . unHelpful . entityDiff $ opts
        vertexContents =
            fmap (fmap (processCsv . CSV.decodeByName) . CL.readFile)
                . unHelpful
                . vertexInput
                $ opts

    hPutStrLn stderr $ "Level information (Name, Number of entities):"
    hPutStrLn stderr
        . show
        . fmap (L.over L._2 (Map.size . unLevel))
        $ levels

    when (isJust vertexContents)
        $ hPutStrLn stderr "Getting vertex similarities."

    vertexSimMap <- maybe
                        (return . defVertexSimMap size $ levelNames)
                        (fmap (vertexCsvToLevels idMap . V.toList))
                  $ vertexContents

    let edgeSimMethod = maybe SpearmanCorrelation read
                      . unHelpful
                      . edgeMethod
                      $ opts
        -- getSimMat ARACNE = getSimMatAracneR size
        -- getSimMat KendallCorrelation = getSimMatKendallR size
        getSimMat SpearmanCorrelation = return . getSimMatCorrelation edgeSimMethod
        getSimMat PearsonCorrelation  = return . getSimMatCorrelation edgeSimMethod
        --getSimMat SpearmanCorrelation = getSimMatSpearmanR size
        -- getSimMat KendallCorrelation = getSimMatKendall
        --                                 eDiff
        --                                 (MaximumEdge 1)
        --                                 idMap

    liftIO $ hPutStrLn stderr "Getting edge similarities."

    edgeSimMap   <- fmap (EdgeSimMap . Map.fromList)
                  . mapM ( L.sequenceOf L._2
                         . L.over L._2 ( getSimMat edgeSimMethod
                                       . standardizeLevel idMap
                                       )
                         )
                  $ levels

    -- let getGrMap KendallCorrelation = getGrKendall eDiff (MaximumEdge 1) idMap
        --getGrMap ARACNE             = undefined

    -- grMap <- if edgeSimMethod == KendallCorrelation
    --             then liftIO
    --                . fmap (GrMap . Map.fromList)
    --                . mapM ( L.sequenceOf L._2
    --                        . L.over L._2 ( getGrMap edgeSimMethod
    --                                        . standardizeLevel idMap
    --                                        )
    --                        )
    --                $ levels
    --             else return . GrMap $ Map.empty
    let grMap = GrMap Map.empty

    return (Nothing, Just unifiedData, idMap, idVec, vertexSimMap, edgeSimMap, grMap)

-- | Get all of the network info that is pre-made for input into the integration method.
getPremadeIntegrationInput
    :: Options
    -> IO (Maybe (Set.Set ID), Maybe UnifiedData, IDMap, IDVec, VertexSimMap, EdgeSimMap, GrMap)
getPremadeIntegrationInput opts = do

    contents <- maybe getContents readFile
              . unHelpful
              . dataInput
              $ opts

    if unHelpful . test $ opts
        then return
           . getPremadeNetworks
           . L.over L._1 Just
           $ (read contents :: ([String], [([String], [(String, String, Double)])]))
        else return
           . getPremadeNetworks
           . (Nothing,)
           $ (read contents :: [([String], [(String, String, Double)])])

-- | Show the accuracy of a test analysis.
showAccuracy :: Set.Set ID -> IDVec -> V.Vector NodeCorrScoresInfo -> T.Text
showAccuracy truthSet idVec nodeCorrScoresInfo =
    T.pack . show . getAccuracy truthSet idVec $ nodeCorrScoresInfo

main :: IO ()
main = do
    opts <- getRecord "integrate, Gregory W. Schwartz\
                      \ Integrate data from multiple sources to find consistent\
                      \ (or inconsistent) entities."

    (truthSet, unifiedData, idMap, idVec, vertexSimMap, edgeSimMap, grMap) <-
        bool (getIntegrationInput opts) (getPremadeIntegrationInput opts)
            . unHelpful
            . premade
            $ opts

    let alignment =
            maybe CosineSimilarity read . unHelpful . alignmentMethod $ opts
        size      = Size . Map.size . unIDMap $ idMap
        nPerm     =
            Permutations . fromMaybe 1000 . unHelpful . permutations $ opts

    hPutStrLn
        stderr
        "Calculating vertex similarities and bootstraps between networks."

    nodeCorrScoresMap <- case alignment of
        CosineSimilarity ->
            integrateCosineSim nPerm size vertexSimMap edgeSimMap
        RandomWalker     ->
            integrateWalker
                nPerm
                size
                ( WalkerRestart
                . fromMaybe 0.25
                . unHelpful
                . walkerRestart
                $ opts
                )
                edgeSimMap
        RandomWalkerSim ->
            integrateWalkerSim
                ( WalkerRestart
                . fromMaybe 0.25
                . unHelpful
                . walkerRestart
                $ opts
                )
                (Counter . fromMaybe 100 . unHelpful . steps $ opts)
                grMap
        -- CSRW -> liftIO
        --     $ integrateCSRW
        --         vertexSimMap
        --         edgeSimMap
        --         ( WalkerRestart
        --         . fromMaybe 0.05
        --         . unHelpful
        --         . walkerRestart
        --         $ opts
        --         )
        --         (Counter . fromMaybe 100 . unHelpful . steps $ opts)

    hPutStrLn stderr "Calculating node correspondence scores."

    nodeCorrScoresInfo <- getNodeCorrScoresInfo nodeCorrScoresMap

    if unHelpful . test $ opts
        then T.putStr
            . maybe
                    (error "Truth set not found.")
                    (\x -> showAccuracy x idVec nodeCorrScoresInfo)
            $ truthSet
        else T.putStr
            . printNodeCorrScores idVec unifiedData nodeCorrScoresMap
            $ nodeCorrScoresInfo

    return ()
