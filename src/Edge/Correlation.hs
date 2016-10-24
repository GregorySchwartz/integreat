{- Edge.Correlation
Gregory W. Schwartz

Collections the functions pertaining to the generation of edges in a similarity
matrix using different kinds of correlations as weights.
-}

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE QuasiQuotes #-}

module Edge.Correlation
    ( getGrKendall
    , getSimMatKendall
    ) where

-- Standard
import Data.Maybe
import Data.Bool
import Data.List
import qualified Data.Sequence as Seq
import qualified Data.Map.Strict as Map
import qualified Data.Foldable as F
import Control.Monad

-- Cabal
import qualified Data.Vector as V
import Data.Graph.Inductive
import qualified Control.Lens as L
import qualified Numeric.LinearAlgebra as N

import qualified Foreign.R as R
import Language.R.Instance as R
import Language.R.QQ
import qualified Language.R.Literal as R

-- Local
import Types
import Utility

-- | Take one level and get the graph by using correlations (a
-- co-expression network). The default value is applied to missing data.
-- Consider using a value that would mean no correlation: 0.
-- Also, possibly correct for added
-- positions in the entity. Use the maximum value if the two entities are the
-- same without the entityDiff.
getGrKendall :: Maybe EntityDiff
             -> MaximumEdge
             -> IDMap
             -> StandardLevel
             -> R s LevelGr
getGrKendall
    entityDiff
    (MaximumEdge maxEdge)
    (IDMap idMap)
    (StandardLevel level) =
    fmap ( LevelGr
         . undir
         . mkGraph (zip [0..Map.size idMap] [0..Map.size idMap])
         . fmap (\((!x, !y), !z) -> (x, y, z))
         . catMaybes
         )
        . pairsM getCorr
        . fmap (L.over L._2 (V.fromList . F.toList))
        . Map.toList
        $ level
  where
      getCorr ((!id1, !idx1), !x) ((!id2, !idx2), !y) = do
          res <- bool (kendallCorrelate x y) (return $ Just maxEdge)
               . sameWithEntityDiff entityDiff id1
               $ id2
          return . fmap ((idx1 , idx2),) $ res
      keyNotFound k = "ID: " ++ show k ++ " not found."

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
                 -> R s EdgeSimMatrix
getSimMatKendall (Default def)
                 entityDiff
                 (MaximumEdge maxEdge)
                 (IDMap idMap)
                 (StandardLevel level) =
    fmap ( EdgeSimMatrix
         . N.assoc (Map.size idMap, Map.size idMap) def
         . mappend diagonal
         . concatMap flipToo
         . catMaybes
         )
        . pairsM getCorr
        . fmap (L.over L._2 (V.fromList . F.toList))
        . Map.toList
        $ level
  where
      diagonal = take (Map.size idMap)
               . flip zip [1,1..]
               . iterate (L.over L.both (+ 1))
               $ (0, 0)
      getCorr ((!id1, !idx1), !x) ((!id2, !idx2), !y) = do
          res <- bool (kendallCorrelate x y) (return $ Just maxEdge)
               . sameWithEntityDiff entityDiff id1
               $ id2
          return . fmap ((idx1 , idx2),) $ res
      keyNotFound k = "ID: " ++ show k ++ " not found."

-- | Correlate two groups of entities, where each group is a collection of
-- measurements (specifically for a single type of entity, for instance
-- a single protein). If there is missing data, we just say their is no
-- correlation: 0.
kendallCorrelate :: V.Vector (Maybe Entity)
                 -> V.Vector (Maybe Entity)
                 -> R s (Maybe Double)
kendallCorrelate e1 e2 = do
    let joined = catMaybes . V.toList . V.zipWith groupDataSets e1 $ e2
        xs     = fmap fst joined
        ys     = fmap snd joined

    if length joined < 2
        then return Nothing
        else do
            res <- [r| cor.test(xs_hs, ys_hs, method="kendall") |]
            p   <- [r| res_hs$p.value |]
            t   <- [r| res_hs$estimate |]

            if (R.fromSomeSEXP p :: Double) >= 0.05
                then return Nothing
                else return . Just $ (R.fromSomeSEXP t :: Double)
