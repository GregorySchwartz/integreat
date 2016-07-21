{- integrate
Gregory W. Schwartz

Integrate data from multiple sources to find consistent (or inconsistent)
entities.
-}

{-# LANGUAGE BangPatterns #-}

module Main where

-- Standard
import qualified Data.Set as Set
import qualified Data.Map.Strict as Map
import Data.Monoid
import System.IO

-- Cabal
import qualified Data.Vector as V
import qualified Data.ByteString.Lazy.Char8 as CL
import qualified Data.Text as T
import qualified Data.Text.IO as T
import qualified Data.Csv as CSV
import Control.Lens
import Options.Applicative

-- Local
import Types
import Utility
import Load
import NetworkGeneration
import RandomWalk
import Integrate
import Print

-- | Command line arguments
data Options = Options { dataInput     :: String
                       , vertexInput   :: Maybe String
                       , method        :: Method
                       , walkerRestart :: Double
                       , counterStop   :: Int
                       }

-- | Command line options
options :: Parser Options
options = Options
      <$> strOption
          ( long "data-input"
         <> short 'd'
         <> metavar "FILE"
         <> help "The input file containing the data intentisities.\
                 \ Follows the format:\
                 \ dataLevel,dataReplicate,vertex,intensity.\
                 \ dataLevel is the level (the base title for the experiment,\
                 \ \"data set\"), dataReplicate is the replicate in that\
                 \ experiment that the entity is\
                 \ from, vertex is the name of the entity (must match those in\
                 \ the vertex-input), and the intensity is the value of this\
                 \ entity in this data set."
          )
      <*> optional ( strOption
          ( long "vertex-input"
         <> short 'v'
         <> metavar "FILE"
         <> help "The input file containing similarities between entities.\
                 \ Follows the format: vertexLevel1,vertexLevel2,\
                 \vertex1,vertex2,similarity.\
                 \ vertexLevel1 is the level (the base title for the\
                 \ experiment, \"data set\") that vertex1 is\
                 \ from, vertexLevel2 is the level that vertex2 is from,\
                 \ and the similarity is a number representing the similarity\
                 \ between those two entities. If not specified, then the same\
                 \ entity (determined by vertex in data-input) will have a\
                 \ similarity of 1, different entities will have a similarity\
                 \ of 0."
          )
        )
      <*> option auto
          ( long "method"
         <> short 'm'
         <> metavar "[CosineSimilarity] | RandomWalker"
         <> help "The method to get integrated vertex similarity between\
                 \ levels. CosineSimilarity uses the cosine similarity of each\
                 \ vertex in each network compared to the other vertices in\
                 \ other networks. RandomWalker uses a random walker based\
                 \ network alignment algorithm in order to get similarity."
          )
      <*> option auto
          ( long "walker-restart"
         <> short 'r'
         <> metavar "[0.05] | PROBABILITY"
         <> value 0.05
         <> help "For the random walker algorithm, the probability of making\
                 \ a jump to a random vertex. Recommended to be the ratio of\
                 \ the total number of vertices in the top 99% smallest\
                 \ subnetworks to the total number of nodes in the reduced\
                 \ product graph (Jeong, 2015)."
          )
      <*> option auto
          ( long "walker-steps"
         <> short 's'
         <> metavar "[10000] | STEPS"
         <> value 10000
         <> help "For the random walker algorithm, the number of steps to take\
                 \ before stopping."
          )

mainFunc :: Options -> IO ()
mainFunc opts = do
    let processCsv = snd . either error id

    dataEntries   <- fmap (\x -> processCsv
                               $ ( CSV.decodeByName x
                                :: Either String (CSV.Header, V.Vector DataEntry)
                                 )
                          )
                   . CL.readFile
                   . dataInput
                   $ opts

    let levels       =
            entitiesToLevels . datasToEntities . V.toList $ dataEntries
        unifiedData  = unifyAllLevels . fmap snd $ levels
        levelNames   = Set.toList . Set.fromList . map fst $ levels
        idMap        = getIDMap unifiedData
        idVec        = getIDVec unifiedData

    vertexSimMap <- maybe
                        (return . defVertexSimMap idMap $ levelNames)
                        (fmap (vertexCsvToLevels idMap . V.toList))
                  . fmap ( fmap ( \x -> processCsv
                                      $ ( CSV.decodeByName x
                                       :: Either
                                            String
                                            (CSV.Header, V.Vector VertexEntry)
                                        )
                                )
                         . CL.readFile
                         )
                  . vertexInput
                  $ opts

    let edgeSimMap = EdgeSimMap
                   . Map.fromList
                   . fmap (over _2 (getSimMat (Default (-5)) idMap))
                   $ levels

    nodeCorrScores <- integrate
                        (method opts)
                        vertexSimMap
                        edgeSimMap
                        (WalkerRestart . walkerRestart $ opts)
                        (Counter . counterStop $ opts)

    T.putStr . printNodeCorrScores idVec $ nodeCorrScores

    return ()

main :: IO ()
main = execParser opts >>= mainFunc
  where
    opts = info (helper <*> options)
      ( fullDesc
     <> progDesc "Integrate data from multiple sources to find consistent\
                 \ (or inconsistent) entities."
     <> header "integrate, Gregory W. Schwartz" )
