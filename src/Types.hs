{- Types
Gregory W. Schwartz

Collections the types used in the program
-}

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE DuplicateRecordFields #-}

module Types where

-- Standard
import qualified Data.Sequence as Seq
import qualified Data.Map.Strict as Map
import qualified Data.IntMap.Strict as IMap
import GHC.Generics

-- Cabal
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Unboxed as VU
import qualified Data.Text as T
import Data.Csv
import Data.Graph.Inductive
import Control.Monad.Reader
import Control.Monad.State
import Control.Lens
import Numeric.LinearAlgebra

-- Local


-- Algebraic
newtype Threads          = Threads Int
newtype ID               = ID { unID :: T.Text } deriving (Eq, Ord, Show)
newtype Default          = Default { unDefault :: Double }
newtype Bandwidth        = Bandwidth { unBandwidth :: Double }
newtype MaximumEdge      = MaximumEdge { unMaximumEdge :: Double }
newtype StdDevThreshold  = StdDevThreshold { unStdDevThreshold :: Double }
newtype EntityDiff       = EntityDiff  { unEntityDiff  :: T.Text }
newtype Size             = Size Int
newtype Counter          = Counter Int deriving (Eq, Ord, Num)
newtype NumSamples       = NumSamples { unNumSamples :: Int }
newtype WalkerRestart    = WalkerRestart { unWalkerRestart :: Double }
newtype DataSetName      = DataSetName T.Text deriving (Eq, Ord, Show)
newtype PValue           = PValue { unPValue :: Double } deriving (Eq, Ord, Show)
newtype Permutations     = Permutations Int
newtype EdgeValues       = EdgeValues (VU.Vector Double)
newtype LevelName        = LevelName { unLevelName :: T.Text }
                           deriving (Eq, Ord, Show)
newtype IDVec            = IDVec { unIDVec :: V.Vector ID }
newtype IDMap            = IDMap { unIDMap :: Map.Map ID Int }
newtype WalkerState      =
    WalkerState { unWalkerState :: (Int, IMap.IntMap Int) }
-- newtype WalkerStateCSRW  =
--     WalkerStateCSRW { unWalkerStateCSRW :: (Int, Int, TransProbMatrix) }

newtype DataSet          = DataSet (Map.Map ID Entity)
newtype StandardDataSets = StandardDataSets
    { unStandardDataSets :: Seq.Seq DataSetName
    }
newtype Level = Level
    { unLevel :: (Map.Map ID (Map.Map DataSetName Entity))
    }
newtype StandardLevel = StandardLevel
    { unStandardLevel :: (Map.Map (ID, Int) (Seq.Seq (Maybe Entity)))
    }
newtype UnifiedData      =
    UnifiedData { unUnifiedData :: (Map.Map ID (Map.Map DataSetName Entity)) }

newtype SimVector =
    SimVector { unSimVector :: VS.Vector Double }
newtype EdgeSimMatrix =
    EdgeSimMatrix { unEdgeSimMatrix :: IMap.IntMap (IMap.IntMap Double) }
    deriving (Show)
newtype VertexSimValues =
    VertexSimValues { unVertexSimValues :: [((Int, Int), Double)] }
    deriving (Show)
newtype TransProbMatrix  =
    TransProbMatrix { unTransProbMatrix :: Matrix Double }
newtype LevelGr = LevelGr { unLevelGr :: Gr Int Double } deriving (Show)
newtype NodeCorrScores   =
    NodeCorrScores { unNodeCorrScores :: V.Vector (Double, Maybe PValue) }
newtype FlatNodeCorrScores =
    FlatNodeCorrScores { unFlatNodeCorrScores :: V.Vector Double }
newtype PValNodeCorrScores =
    PValNodeCorrScores { unPValNodeCorrScores :: V.Vector Double }
newtype NodeCorrScoresMap = NodeCorrScoresMap
    { unNodeCorrScoresMap :: Map.Map (LevelName, LevelName) NodeCorrScores
    }

newtype EdgeSimMap       =
    EdgeSimMap { unEdgeSimMap :: (Map.Map LevelName EdgeSimMatrix) }
    deriving (Show)
newtype GrMap            =
    GrMap { unGrMap :: (Map.Map LevelName LevelGr) }
    deriving (Show)
newtype VertexSimMap     =
    VertexSimMap { unVertexSimMap
                :: (Map.Map LevelName (Map.Map LevelName VertexSimValues))
                 }
    deriving (Show)

newtype Walker a =
    Walker { unWalker :: (ReaderT Environment (StateT WalkerState IO) a) }
    deriving ( Functor
             , Applicative
             , Monad
             , MonadIO
             , MonadReader Environment
             , MonadState WalkerState
             )
-- newtype WalkerCSRW a =
--     WalkerCSRW { unWalkerCSRW :: (ReaderT EnvironmentCSRW (StateT WalkerStateCSRW IO) a) }
--     deriving ( Functor
--              , Applicative
--              , Monad
--              , MonadIO
--              , MonadReader EnvironmentCSRW
--              , MonadState WalkerStateCSRW
--              )

data Entity = Entity { _entityID    :: !ID
                     , _dataSetName :: !DataSetName
                     , _levelName   :: !LevelName
                     , _entityValue :: !Double
                     }
              deriving (Eq, Ord, Show)

data Environment =
    Environment { eGr     :: !LevelGr
                , restart :: !WalkerRestart
                , v0      :: !Int
                }
-- data EnvironmentCSRW =
--     EnvironmentCSRW { edgeMat1  :: !EdgeSimMatrix
--                     , edgeMat2  :: !EdgeSimMatrix
--                     , vertexMat :: !VertexSimMatrix
--                     , restart   :: !WalkerRestart
--                     }

data WalkerChoice    = Same | DifferentLeft | DifferentRight
data AlignmentMethod
    = CosineSimilarity
    | RandomWalker
    deriving (Eq,Read,Show)
    -- | CSRW
data EdgeMethod = ARACNE | SpearmanCorrelation | KendallCorrelation deriving (Eq,Read,Show)

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

data NodeCorrScoresInfo = NodeCorrScoresInfo
    { nodeCorrScoresMap          :: NodeCorrScoresMap
    , avgNodeCorrScores          :: FlatNodeCorrScores
    , avgPValNodeCorrScores      :: PValNodeCorrScores
    , rankProdNodeCorrScores     :: FlatNodeCorrScores
    , rankProdPValNodeCorrScores :: PValNodeCorrScores
    }

makeLenses ''Entity

-- Basic

-- Advanced
