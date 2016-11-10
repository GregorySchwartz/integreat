{- Alignment.Cosine
Gregory W. Schwartz

Collections the functions pertaining to getting integration by cosine
similarity.
-}

{-# LANGUAGE BangPatterns #-}

module Alignment.Cosine
    ( cosineIntegrate
    ) where

-- Standard
import Data.Tuple
import Data.List
import qualified Data.IntMap.Strict as IMap

-- Cabal
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import Control.Lens
import Numeric.LinearAlgebra

-- Local
import Types
import Utility

-- | Get the cosine similarity of two levels. Normalize based on the size
-- of actual data (no missing data) divided by the total amount of data.
cosineIntegrate :: Size
                -> VertexSimMap
                -> LevelName
                -> LevelName
                -> EdgeSimMatrix
                -> EdgeSimMatrix
                -> NodeCorrScores
cosineIntegrate size vMap l1 l2 e1 e2 =
    NodeCorrScores
        . V.fromList
        . fmap snd
        . IMap.toAscList
        . fillIntMap size
        . IMap.intersectionWith cosineSimIMap newE1
        $ newE2
  where
    newE1 :: IMap.IntMap (IMap.IntMap Double)
    newE1         =
        unEdgeSimMatrix . cosineUpdateSimMat e1 . unVertexSimValues $ vertexSim
    newE2         =
        unEdgeSimMatrix . cosineUpdateSimMat e2 . unVertexSimValues $ vertexSim
    vertexSim     = getVertexSim l1 l2 vMap

-- | Update the similarity matrices where the diagonal contains the similarity
-- between vertices of different levels.
cosineUpdateSimMat :: EdgeSimMatrix -> [((Int, Int), Double)] -> EdgeSimMatrix
cosineUpdateSimMat (EdgeSimMatrix edgeSimMat) =
    EdgeSimMatrix
        . foldl' (\acc ((!x, !y), !z) -> IMap.alter (f y z) x acc) edgeSimMat
        . (\xs -> xs ++ fmap (over _1 swap) xs)
  where
    f y z Nothing  = Just . IMap.singleton y $ z
    f y z (Just v) = Just . IMap.insert y z $ v
