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
import Data.Bool
import Data.Maybe
import qualified Data.Set as Set
import qualified Data.Map.Strict as Map
import Data.Monoid
import System.IO
import Control.Monad

-- Cabal
import qualified Data.Vector as V
import qualified Data.ByteString.Lazy.Char8 as CL
import qualified Data.Text as T
import qualified Data.Text.IO as T
import qualified Data.Csv as CSV
import qualified Control.Lens as L
import Control.Monad.Trans
import Options.Generic

import qualified Foreign.R as R
import Language.R.Instance as R
import Language.R.QQ

-- Local
import Types
import Utility
import Load
import Edge.Correlation
import Edge.Aracne
import Alignment.CSRW
import Integrate
import Print

-- | Command line arguments
data Options = Options { dataInput       :: Maybe String
                                        <?> "(FILE) The input file containing the data intensities. Follows the format: dataLevel,dataReplicate,vertex,intensity. dataLevel is the level (the base level for the experiment, like \"proteomic_Myla\" or \"RNA_MyLa\" for instance), dataReplicate is the replicate in that experiment that the entity is from (the name of that data set with the replicate name, like \"RNA_MyLa_1\"), and vertex is the name of the entity (must match those in the vertex-input), and the intensity is the value of this entity in this data set."
                       , vertexInput     :: Maybe String
                                        <?> "(FILE) The input file containing similarities between entities. Follows the format: vertexLevel1,vertexLevel2, vertex1,vertex2,similarity. vertexLevel1 is the level (the base title for the experiment, \"data set\") that vertex1 is from, vertexLevel2 is the level that vertex2 is from, and the similarity is a number representing the similarity between those two entities. If not specified, then the same entity (determined by vertex in data-input) will have a similarity of 1, different entities will have a similarity of 0."
                       , entityDiff      :: Maybe T.Text
                                        <?> "When comparing entities that are the same, ignore the text after this separator. Used for comparing phosphorylated positions with another level. For example, if we have a strings ARG29 and ARG29_7 that we want to compare, we want to say that their value is the highest in correlation, so this string would be \"_\""
                       , alignmentMethod :: Maybe String
                                        <?> "([CosineSimilarity] | RandomWalker | CSRW) The method to get integrated vertex similarity between levels. CosineSimilarity uses the cosine similarity of each  vertex in each network compared to the other vertices in  other networks. RandomWalker uses a random walker based  network alignment algorithm in order to get similarity."
                       , edgeMethod      :: Maybe String
                                        <?> "([KendallCorrelation] | ARACNE) The method to use for the edges between entities in the coexpression matrix."
                       , walkerRestart   :: Maybe Double
                                        <?> "([0.25] | PROBABILITY) For the random walker algorithm, the probability of making  a jump to a random vertex. Recommended to be the ratio of  the total number of vertices in the top 99% smallest  subnetworks to the total number of nodes in the reduced  product graph (Jeong, 2015)."
                       , steps           :: Maybe Int
                                        <?> "([100] | STEPS) For the random walker algorithm, the number of steps to take  before stopping."
                       , premade         :: Bool
                                        <?> "Whether the input data (dataInput) is a pre-made network of the format \"[([\"VERTEX\"], [(\"SOURCE\", \"DESTINATION\", WEIGHT)])]\", where VERTEX, SOURCE, and DESTINATION are of type INT starting at 0, in order, and WEIGHT is a DOUBLE representing the weight of the edge between SOURCE and DESTINATION."
                       , test            :: Bool
                                        <?> "Whether the input data from premade is from a test run. If supplied, the output is changed to an accuracy measure. In this case, we get the total rank below the number of permuted vertices divided by the theoretical maximum (so if there were five changed vertices out off 10 and two were rank 8 and 10 while the others were in the top five, we would have (1 - ((3 + 5) / (10 + 9 + 8 + 7 + 6))) as the accuracy."
                       , entityFilter    :: Maybe Int
                                        <?> "The minimum number of samples an entity must appear in, otherwise the entity is removed from the analysis."
                       }
               deriving (Generic)

instance ParseRecord Options

-- | Get all of the required information for integration.
getIntegrationInput
    :: Options
    -> R s (Maybe (Set.Set ID), Maybe UnifiedData, IDMap, IDVec, VertexSimMap, EdgeSimMap, GrMap)
getIntegrationInput opts = do
    let processCsv = snd . either error id

    liftIO $ hPutStrLn stderr "Getting data input."
    dataEntries   <- liftIO
                   . fmap (processCsv . CSV.decodeByName)
                   . maybe CL.getContents CL.readFile
                   . unHelpful
                   . dataInput
                   $ opts

    let levels         = entitiesToLevels
                       . (\x -> maybe x (flip filterEntities x) numSamples)
                       . datasToEntities
                       . V.toList
                       $ dataEntries
        numSamples     = fmap NumSamples . unHelpful . entityFilter $ opts
        unifiedData    = unifyAllLevels . fmap snd $ levels
        levelNames     = Set.toList . Set.fromList . fmap fst $ levels
        idMap          = getIDMap unifiedData
        idVec          = getIDVec unifiedData
        eDiff          = fmap EntityDiff . unHelpful . entityDiff $ opts
        vertexContents =
            fmap (fmap (processCsv . CSV.decodeByName) . CL.readFile)
                . unHelpful
                . vertexInput
                $ opts

    when (isJust vertexContents)
        . liftIO
        $ hPutStrLn stderr "Getting vertex similarities."
        
    vertexSimMap <- liftIO
                  . maybe
                        (return . defVertexSimMap idMap $ levelNames)
                        (fmap (vertexCsvToLevels idMap . V.toList))
                  $ vertexContents

    let edgeSimMethod = maybe KendallCorrelation read
                      . unHelpful
                      . edgeMethod
                      $ opts
        getSimMat ARACNE = getSimMatAracneR
        getSimMat KendallCorrelation = getSimMatKendallR
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

    let getGrMap KendallCorrelation = getGrKendall eDiff (MaximumEdge 1) idMap
        getGrMap ARACNE             = undefined

    grMap <- if edgeSimMethod == KendallCorrelation
                then liftIO
                   . fmap (GrMap . Map.fromList)
                   . mapM ( L.sequenceOf L._2
                           . L.over L._2 ( getGrMap edgeSimMethod
                                           . standardizeLevel idMap
                                           )
                           )
                   $ levels
                else return . GrMap $ Map.empty

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
showAccuracy :: Set.Set ID -> IDVec -> NodeCorrScoresInfo -> T.Text
showAccuracy truthSet idVec nodeCorrScoresInfo =
    T.pack . show . getAccuracy truthSet idVec $ nodeCorrScoresInfo

main :: IO ()
main = do
    opts <- getRecord "integrate, Gregory W. Schwartz\
                      \ Integrate data from multiple sources to find consistent\
                      \ (or inconsistent) entities."

    hPutStrLn stderr "Starting R thread within Haskell."
    
    R.withEmbeddedR R.defaultConfig $ R.runRegion $ do
        (truthSet, unifiedData, idMap, idVec, vertexSimMap, edgeSimMap, grMap) <-
            bool (getIntegrationInput opts) (liftIO $ getPremadeIntegrationInput opts)
                . unHelpful
                . premade
                $ opts

        let alignment =
                maybe CosineSimilarity read . unHelpful . alignmentMethod $ opts

        liftIO $
            hPutStrLn stderr "Calculating vertex similarities between networks."
        
        nodeCorrScoresMap <- case alignment of
            CosineSimilarity -> liftIO
                $ integrateCosineSim vertexSimMap edgeSimMap
            RandomWalker -> liftIO
                $ integrateWalker
                    ( WalkerRestart
                    . fromMaybe 0.25
                    . unHelpful
                    . walkerRestart
                    $ opts
                    )
                    (Counter . fromMaybe 100 . unHelpful . steps $ opts)
                    grMap
            CSRW -> liftIO
                $ integrateCSRW
                    vertexSimMap
                    edgeSimMap
                    ( WalkerRestart
                    . fromMaybe 0.05
                    . unHelpful
                    . walkerRestart
                    $ opts
                    )
                    (Counter . fromMaybe 100 . unHelpful . steps $ opts)

        liftIO $ hPutStrLn stderr "Calculating node correspondence scores."
            
        nodeCorrScoresInfo <- getNodeCorrScoresInfo nodeCorrScoresMap

        if unHelpful . test $ opts
            then liftIO
               . T.putStr
               . maybe
                       (error "Truth set not found.")
                       (\x -> showAccuracy x idVec nodeCorrScoresInfo)
               $ truthSet
            else liftIO
               . T.putStr
               . printNodeCorrScores idVec unifiedData
               $ nodeCorrScoresInfo

        return ()
