{- Types
Gregory W. Schwartz

Collections the types used in the program
-}

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE TemplateHaskell #-}

module Types where

-- Standard
import qualified Data.Map.Strict as Map
import qualified Data.Sequence as Seq
import GHC.Generics

-- Cabal
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Data.Text as T
import Data.Csv
import Control.Monad.Reader
import Control.Monad.State
import Control.Monad.Random
import Control.Lens
import Numeric.LinearAlgebra

-- Local


-- Algebraic
newtype ID               = ID { unID :: T.Text } deriving (Eq, Ord, Show)
newtype Default          = Default { unDefault :: Double }
newtype Size             = Size Int
newtype Counter          = Counter Int deriving (Eq, Ord, Num)
newtype WalkerRestart    = WalkerRestart { unWalkerRestart :: Double }
newtype DataSetName      = DataSetName T.Text deriving (Show)
newtype LevelName        = LevelName { unLevelName :: T.Text }
                           deriving (Eq, Ord, Show)
newtype IDVec            = IDVec { unIDVec :: V.Vector ID }
newtype IDMap            = IDMap { unIDMap :: Map.Map ID Int }
newtype WalkerState      =
    WalkerState { unWalkerState :: (Int, Int, TransProbMatrix) }

newtype DataSet          = DataSet (Map.Map ID Entity)
newtype Level            = Level { unLevel :: (Map.Map ID (Seq.Seq Entity)) }
newtype UnifiedData      =
    UnifiedData { unUnifiedData :: (Map.Map ID (Seq.Seq Entity)) }

newtype SimVector =
    SimVector { unSimVector :: VS.Vector Double }
newtype EdgeSimMatrix =
    EdgeSimMatrix { unEdgeSimMatrix :: Matrix Double }
    deriving (Show)
newtype VertexSimMatrix =
    VertexSimMatrix { unVertexSimMatrix :: Matrix Double }
    deriving (Show)
newtype TransProbMatrix  =
    TransProbMatrix { unTransProbMatrix :: Matrix Double }
newtype NodeCorrScores   =
    NodeCorrScores { unNodeCorrScores :: VS.Vector Double }

newtype EdgeSimMap       =
    EdgeSimMap { unEdgeSimMap :: (Map.Map LevelName EdgeSimMatrix) }
    deriving (Show)
newtype VertexSimMap     =
    VertexSimMap { unVertexSimMap
                :: (Map.Map LevelName (Map.Map LevelName VertexSimMatrix))
                 }
    deriving (Show)

newtype Walker a =
    Walker { unWalker :: (ReaderT Environment (StateT WalkerState IO) a) }
    deriving ( Functor
             , Applicative
             , Monad
             , MonadIO
             , MonadRandom
             , MonadReader Environment
             , MonadState WalkerState
             )

data Entity = Entity { _entityID    :: !ID
                     , _dataSetName :: !DataSetName
                     , _levelName   :: !LevelName
                     , _entityValue :: !Double
                     }
              deriving (Show)

data Environment =
    Environment { edgeMat1  :: !EdgeSimMatrix
                , edgeMat2  :: !EdgeSimMatrix
                , vertexMat :: !VertexSimMatrix
                , restart   :: !WalkerRestart
                }

data WalkerChoice = Same | DifferentLeft | DifferentRight
data Method       = CosineSimilarity | RandomWalker deriving (Eq, Read, Show)

data DataEntry    = DataEntry { dataLevel     :: !T.Text
                              , dataReplicate :: !T.Text
                              , vertex        :: !T.Text
                              , intensity     :: !Double
                              }
                    deriving (Generic)

instance FromNamedRecord DataEntry
instance ToNamedRecord DataEntry
instance DefaultOrdered DataEntry

data VertexEntry  = VertexEntry { vertexLevel1 :: !T.Text
                                , vertexLevel2 :: !T.Text
                                , vertex1      :: !T.Text
                                , vertex2      :: !T.Text
                                , similarity   :: !Double
                                }
                    deriving (Generic)

instance FromNamedRecord VertexEntry
instance ToNamedRecord VertexEntry
instance DefaultOrdered VertexEntry

makeLenses ''Entity

-- Basic

-- Advanced
