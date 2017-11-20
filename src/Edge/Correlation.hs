{- Edge.Correlation
Gregory W. Schwartz

Collections the functions pertaining to the generation of edges in a similarity
matrix using different kinds of correlations as weights.
-}

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE GADTs #-}

module Edge.Correlation
    ( getSimMatCorrelation
    , correlate
    ) where

-- Standard
import Control.Monad
import Data.Bool
import Data.List
import Data.Maybe
import qualified Data.Foldable as F
import qualified Data.IntMap.Strict as IMap
import qualified Data.Map.Strict as Map
import qualified Data.Sequence as Seq
import qualified Data.Set as Set

-- Cabal
import Control.Monad.Trans
import Data.Graph.Inductive
import Statistics.Correlation
import Statistics.Distribution
import Statistics.Distribution.StudentT
import qualified Control.Foldl as Fold
import qualified Control.Lens as L
import qualified Data.Vector.Unboxed as VU
import qualified Numeric.LinearAlgebra as N

-- Local
import Types
import Utility

getSimMatCorrelation :: EdgeMethod -> StandardLevel -> EdgeSimMatrix
getSimMatCorrelation edgeMethod (StandardLevel level) = do
    EdgeSimMatrix
        . Fold.fold (Fold.Fold step begin extract)
        . pairs (,)
        . fmap (L.over L._2 F.toList)
        . Map.toList
        $ level
  where
      n = fromIntegral . Map.size $ level
      step acc x =
        case uncurry getCorr x of
            ((_, _), Nothing) -> acc
            ((!idx1, !idx2), (Just !val)) ->
                alterMap idx1 idx2 val
                    . alterMap idx2 idx1 val
                    $ acc
      begin      = IMap.empty
      extract    = id
      alterMap new query val =
            IMap.alter ( maybe (Just . IMap.singleton new $ val)
                               (Just . IMap.insert new val)
                       )
                        query
      getCorr ((!id1, !idx1), !x) ((!id2, !idx2), !y) =
            ((idx1, idx2), correlateEntities edgeMethod x y)

-- | Correlate two groups of entities, where each group is a collection of
-- measurements (specifically for a single type of entity, for instance
-- a single protein).
correlateEntities :: EdgeMethod
                  -> [Maybe Entity]
                  -> [Maybe Entity]
                  -> Maybe Double
correlateEntities edgeMethod e1 e2 = do
    let joined = VU.fromList . catMaybes . zipWith groupDataSets e1 $ e2

    guard $ (VU.length joined > 2)
         && (not . VU.all (== (fst . VU.head $ joined)) . VU.map fst $ joined)
         && (not . VU.all (== (snd . VU.head $ joined)) . VU.map snd $ joined)

    let (Rho !coeff, P !pVal) = correlate edgeMethod joined

    guard $ pVal < 0.05

    return coeff

-- | Correlate two lists with p values.
correlate :: EdgeMethod -> VU.Vector (Double, Double) -> (Rho, P)
correlate method xs =
    (Rho (coeff method), P (pVal s))
  where
    n     = fromIntegral . VU.length $ xs
    coeff SpearmanCorrelation = spearman xs
    coeff PearsonCorrelation  = pearson xs
    s     = coeff method * ((sqrt (n - 2)) / (1 - (coeff method ^ 2)))
    pVal stat = 2 * (complCumulative (studentT (fromIntegral $ VU.length xs - 2)) . abs $ stat)
