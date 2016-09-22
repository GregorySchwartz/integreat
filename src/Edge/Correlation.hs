{- Edge.Correlation
Gregory W. Schwartz

Collections the functions pertaining to the generation of edges in a similarity
matrix using different kinds of correlations as weights.
-}

{-# LANGUAGE BangPatterns #-}

module Edge.Correlation
    ( getSimMatKendall
    ) where

-- Standard
import Data.Maybe
import Data.Bool
import Data.List
import qualified Data.Sequence as Seq
import qualified Data.Map.Strict as Map
import qualified Data.Foldable as F

-- Cabal
import qualified Data.Vector as V
import qualified Control.Lens as L
import Numeric.LinearAlgebra
import Statistics.Correlation.Kendall

-- Local
import Types
import Utility

-- | Take one level and get the similarity matrix by using correlations (a
-- co-expression network). The default value is applied to missing data.
-- Consider using a value that would mean no correlation: 0.
-- Also, possibly correct for added
-- positions in the entity. Use the maximum value if the two entities are the
-- same without the entityDiff.
getSimMatKendall :: Default
                 -> Maybe EntityDiff
                 -> MaximumEdge
                 -> IDMap
                 -> StandardLevel
                 -> EdgeSimMatrix
getSimMatKendall (Default def)
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
          , bool (kendallCorrelate x y) maxEdge
          . sameWithEntityDiff entityDiff id1
          $ id2
          )
      keyNotFound k = "ID: " ++ show k ++ " not found."

-- | Correlate two groups of entities, where each group is a collection of
-- measurements (specifically for a single type of entity, for instance
-- a single protein). If there is missing data, we just say their is no
-- correlation: 0.
kendallCorrelate :: V.Vector (Maybe Entity) -> V.Vector (Maybe Entity) -> Double
kendallCorrelate e1 =
    (\x -> if isNaN x then 0 else x)
        . kendall
        . V.fromList
        . catMaybes
        . V.toList
        . V.zipWith groupDataSets e1
