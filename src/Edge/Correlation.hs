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
    -- , getGrKendall
    -- , getSimMatKendall
    -- , getSimMatKendallR
    -- , getSimMatSpearmanR
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

-- import qualified Foreign.R as R
-- import Language.R.Instance as R
-- import Language.R.QQ
-- import qualified Language.R.Literal as R

-- Local
import Types
import Utility

-- | Take one level and get the graph by using correlations (a
-- co-expression network). The default value is applied to missing data.
-- Consider using a value that would mean no correlation: 0.
-- Also, possibly correct for added
-- positions in the entity. Use the maximum value if the two entities are the
-- same without the entityDiff.
-- getGrKendall :: Maybe EntityDiff
--              -> MaximumEdge
--              -> IDMap
--              -> StandardLevel
--              -> IO LevelGr
-- getGrKendall
--     entityDiff
--     (MaximumEdge maxEdge)
--     (IDMap idMap)
--     (StandardLevel level) =
--     fmap (LevelGr . undir)
--         . Fold.foldM (Fold.FoldM step begin extract)
--         . pairs (,)
--         . fmap (L.over L._2 F.toList)
--         . Map.toList
--         $ level
--   where
--     step acc x = do
--       res <- uncurry getCorr x
--       case res of
--           Nothing     -> return acc
--           (Just res') -> return . insEdge res' $ acc
--     begin      =
--       return $ mkGraph (zip [0..Map.size idMap] [0..Map.size idMap]) []
--     extract    = return
--     getCorr ((!id1, !idx1), !x) ((!id2, !idx2), !y) = do
--         res <- bool (kendallCorrelate x y) (return $ Just maxEdge)
--              . sameWithEntityDiff entityDiff id1
--              $ id2
--         return . fmap (idx1, idx2,) $ res
--     keyNotFound k = "ID: " ++ show k ++ " not found."

-- | Take one level and get the similarity matrix by using correlations (a
-- co-expression network). The default value is applied to missing data.
-- Consider using a value that would mean no correlation: 0.
-- Also, possibly correct for added
-- positions in the entity. Use the maximum value if the two entities are the
-- same without the entityDiff.
-- getSimMatKendall :: Maybe EntityDiff
--                  -> MaximumEdge
--                  -> IDMap
--                  -> StandardLevel
--                  -> IO EdgeSimMatrix
-- getSimMatKendall entityDiff
--                  (MaximumEdge maxEdge)
--                  (IDMap idMap)
--                  (StandardLevel level) =
--     fmap EdgeSimMatrix
--         . Fold.foldM (Fold.FoldM step begin extract)
--         . pairs (,)
--         . fmap (L.over L._2 F.toList)
--         . Map.toList
--         $ level
--   where
--       step acc x = do
--         res <- uncurry getCorr x
--         case res of
--             Nothing     -> return acc
--             (Just ((!idx1, !idx2), !val)) ->
--                 return
--                     . alterMap idx1 idx2 val
--                     . alterMap idx2 idx1 val
--                     $ acc
--       begin      = return IMap.empty
--       extract    = return
--       alterMap new query val =
--             IMap.alter ( maybe (Just . IMap.singleton new $ val)
--                                (Just . IMap.insert new val)
--                        )
--                         query
--       getCorr ((!id1, !idx1), !x) ((!id2, !idx2), !y) = do
--           res <- bool (kendallCorrelate x y) (return $ Just maxEdge)
--                . sameWithEntityDiff entityDiff id1
--                $ id2
--           return . fmap (\ !x -> ((idx1, idx2), x)) $ res
--       keyNotFound k = "ID: " ++ show k ++ " not found."

-- | Complete R version.
-- Take one level and get the similarity matrix by using correlations (a
-- co-expression network). The default value is applied to missing data.
-- Consider using a value that would mean no correlation: 0.
-- Does *not* correct for positions with entityDiff here.
-- getSimMatKendallR :: Size -> StandardLevel -> R s EdgeSimMatrix
-- getSimMatKendallR size level = do
--     rDF <- standardLevelToRJSON level

--     [r| suppressPackageStartupMessages(library("psych")) |]
--     rMat <- [r| write("Getting Kendall correlations.", stderr());
--                 df = corr.test(rDF_hs, method = "kendall", adjust = "none", ci = FALSE);
--                 write("Setting bad p-value correlations to 0", stderr());
--                 df$r[df$p >= 0.05] = 0;
--                 df$r
--             |]

--     res <- rToMat size rMat
--     return . EdgeSimMatrix $ res

-- | Correlate two groups of entities, where each group is a collection of
-- measurements (specifically for a single type of entity, for instance
-- a single protein). If there is missing data, we just say there is no
-- correlation: 0.
-- kendallCorrelate :: [Maybe Entity] -> [Maybe Entity] -> IO (Maybe Double)
-- kendallCorrelate e1 e2 = R.runRegion $ do
--     let joined = catMaybes . zipWith groupDataSets e1 $ e2
--         xs     = fmap fst joined
--         ys     = fmap snd joined

--     if length joined < 2
--         then return Nothing
--         else do
--             res <- [r| cor.test(xs_hs, ys_hs, method="kendall") |]
--             p   <- [r| res_hs$p.value |]
--             t   <- [r| res_hs$estimate |]

--             if (R.fromSomeSEXP p :: Double) >= 0.05
--                 then return Nothing
--                 else return . Just $ (R.fromSomeSEXP t :: Double)

-- | Take one level and get the similarity matrix by using correlations (a
-- co-expression network). The default value is applied to missing data.
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

-- | Complete R version.
-- Take one level and get the similarity matrix by using correlations (a
-- co-expression network). The default value is applied to missing data.
-- Consider using a value that would mean no correlation: 0.
-- Does *not* correct for positions with entityDiff here.
-- getSimMatSpearmanR :: Size -> StandardLevel -> R s EdgeSimMatrix
-- getSimMatSpearmanR size level = do
--     rDF <- standardLevelToRJSON level

--     [r| suppressPackageStartupMessages(library("psych")) |]
--     rMat <- [r| write("Getting Spearman correlations.", stderr());
--                 df = corr.test(rDF_hs, method = "spearman", adjust = "none", ci = FALSE);
--                 write("Setting bad p-value correlations to 0", stderr());
--                 df$r[df$p >= 0.05] = 0
--                 df$r
--             |]

--     res <- rToMat size rMat
--     return . EdgeSimMatrix $ res
