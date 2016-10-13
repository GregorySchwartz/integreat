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
cosineIntegrate :: VertexSimMap
                -> LevelName
                -> LevelName
                -> EdgeSimMatrix
                -> EdgeSimMatrix
                -> NodeCorrScores
cosineIntegrate vMap l1 l2 e1 e2 =
    NodeCorrScores
        . VS.convert
        . VS.map applyCosine
        . VS.enumFromN 0
        . rows
        . unEdgeSimMatrix
        $ e1
  where
    applyCosine x = cosineSim (newE1 ! x) (newE2 ! x)
    newE1         = unEdgeSimMatrix . cosineUpdateSimMat e1 $ changes
    newE2         = unEdgeSimMatrix . cosineUpdateSimMat e2 $ changes
    changes       = V.toList
                  . V.imap (\i v -> ((i, i), v))
                  . VS.convert
                  . takeDiag
                  . unVertexSimMatrix
                  $ vertexSim
    vertexSim     = getVertexSim l1 l2 vMap

-- | Update the similarity matrices where the diagonal contains the similarity
-- between vertices of different levels.
cosineUpdateSimMat :: EdgeSimMatrix -> [((Int, Int), Double)] -> EdgeSimMatrix
cosineUpdateSimMat (EdgeSimMatrix edgeSimMat) =
    EdgeSimMatrix
        . accum edgeSimMat const
        . (\xs -> xs ++ fmap (over _1 swap) xs)
