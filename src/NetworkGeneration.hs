{- NetworkGeneration
Gregory W. Schwartz

Collections the functions pertaining to the generation of vertices and
edges in the data sets using a similarity matrix.
-}

{-# LANGUAGE BangPatterns #-}

module NetworkGeneration
    ( getSimMat
    ) where

-- Standard
import Data.Bool
import qualified Data.Map.Strict as Map
import qualified Data.Sequence as Seq
import qualified Data.Foldable as F

-- Cabal
import qualified Data.Vector as V
import qualified Control.Lens as L
import Statistics.Correlation.Kendall
import Numeric.LinearAlgebra

-- Local
import Types
import Utility

-- | Take one level and get the similarity matrix by using correlations (a
-- co-expression network). The default value is applied to missing data.
-- Consider using a value that would never occur: for instance with the
-- Kendall correlation, -5 would never appear. Also, possibly correct for added
-- positions in the entity. Use the maximum value if the two entities are the
-- same without the entityDiff.
getSimMat :: Default
          -> Maybe EntityDiff
          -> MaximumEdge
          -> IDMap
          -> StandardLevel
          -> EdgeSimMatrix
getSimMat (Default def)
          entityDiff
          (MaximumEdge maxEdge)
          (IDMap idMap)
          (StandardLevel level) =
    EdgeSimMatrix
        . assoc (Map.size idMap, Map.size idMap) def
        . concatMap flipToo
        . pairs getCorr
        . fmap (L.over L._2 (V.fromList . F.toList))
        . Map.toList
        $ level
  where
      getCorr ((!id1, !idx1), !x) ((!id2, !idx2), !y) =
          ( (idx1 , idx2)
          , bool (correlate x y) maxEdge
          . sameWithEntityDiff entityDiff id1
          $ id2
          )
      keyNotFound k = "ID: " ++ show k ++ " not found."

-- | Correlate two groups of entities, where each group is a collection of
-- measurements (specifically for a single type of entity, for instance
-- a single protein). If there is missing data, we just say the default
-- value that could never exist: -5.
correlate :: V.Vector Entity -> V.Vector Entity -> Double
correlate e1 e2 =
    (\x -> if isNaN x then (-5) else x)
        . kendall
        . V.zipWith (\ !x !y -> (_entityValue x, _entityValue y)) e1
        $ e2
