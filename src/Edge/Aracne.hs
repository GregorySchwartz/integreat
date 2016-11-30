{- Edge.Aracne
Gregory W. Schwartz

Collections the functions pertaining to the generation of edges using the ARACNE algorithm.

Margolin, A. A., Nemenman, I., Basso, K., Wiggins, C., Stolovitzky, G.,
Favera, R. D., & Califano, A. (). ARACNE: an algorithm for the reconstruction
of gene regulatory networks in a mammalian cellular context. , 7(1), 1–15.
http://dx.doi.org/10.1186/1471-2105-7-S1-S7
-}

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE QuasiQuotes #-}

module Edge.Aracne
    ( getSimMatAracneR
    -- , getSimMatAracne
    ) where

-- Standard
import Data.Maybe
import Data.Bool
import Data.List
import qualified Data.Sequence as Seq
import qualified Data.Map.Strict as Map
import qualified Data.Foldable as F
import Data.Function (on)
import Debug.Trace

-- Cabal
import qualified Data.Vector as V
import qualified Control.Lens as L
import Numeric.LinearAlgebra
import Statistics.Sample

import qualified Foreign.R as R
import Language.R.Instance as R
import Language.R.QQ
import qualified Language.R.Literal as R

-- Local
import Types
import Utility

-- -- | Take one level and get the similarity matrix by using the ARACNE algorithm.
-- -- Also, possibly correct for added positions in the entity. Use the maximum
-- -- value if the two entities are the same without the entityDiff.
-- getSimMatAracne :: Bandwidth
--                 -> Maybe EntityDiff
--                 -> IDMap
--                 -> StandardLevel
--                 -> EdgeSimMatrix
-- getSimMatAracne h entityDiff (IDMap idMap) (StandardLevel level) =
--     dpiPrune
--         . EdgeSimMatrix
--         . assoc (Map.size idMap, Map.size idMap) 0
--         . concatMap flipToo
--         . rankEdges
--         . (\x -> traceShow x x)
--         . filter ((/= 0) . snd)
--         . pairs getMI
--         . fmap (L.over L._2 (V.fromList . F.toList))
--         . Map.toList
--         $ level
--   where
--       getMI ((!id1, !idx1), !x) ((!id2, !idx2), !y) =
--           ( (idx1 , idx2)
--           , bool (aracne h x y) 1
--           . sameWithEntityDiff entityDiff id1
--           $ id2
--           )
--       keyNotFound k = "ID: " ++ show k ++ " not found."

-- -- | Prune the similarity matrix using Data Processing Inequality.
-- dpiPrune :: EdgeSimMatrix -> EdgeSimMatrix
-- dpiPrune (EdgeSimMatrix mat) =
--     EdgeSimMatrix
--         . cmap (\x -> if x /= 0 then 1 else 0)
--         . accum mat const
--         . fmap (, 0)
--         . triples (dpi (EdgeSimMatrix mat))
--         $ [0..(rows mat - 1)]

-- -- | Data Processing Inequality: return the edge to remove if it's smaller than
-- -- the other two.
-- dpi :: EdgeSimMatrix -> Int -> Int -> Int -> (Int, Int)
-- dpi (EdgeSimMatrix mat) g1 g2 g3 =
--     minimumBy (compare `on` (atIndex mat)) [(g1, g3), (g1, g2), (g2, g3)]

-- -- | Rank the edges of a similarity association list.
-- rankEdges :: [((Int, Int), Double)] -> [((Int, Int), Double)]
-- rankEdges = concat
--           . zipWith (\x -> fmap (L.set L._2 x)) [1..]
--           . groupBy ((==) `on` snd)
--           . sortBy (compare `on` snd)

-- -- | Estimate the MI from two random variables, here groups of values from
-- -- samples.
-- aracne :: Bandwidth
--        -> V.Vector (Maybe Entity)
--        -> V.Vector (Maybe Entity)
--        -> Double
-- aracne h e1 = miEstimate h
--             . V.fromList
--             . catMaybes
--             . V.toList
--             . V.zipWith groupDataSets e1

-- -- | The MI estimate.
-- miEstimate :: Bandwidth -> V.Vector (Double, Double) -> Double
-- miEstimate h samples = (1 / (fromIntegral . V.length $ samples))
--                      * (V.sum . V.map (uncurry processSamples) $ samples)
--   where
--     processSamples !x !y = log
--                          $ bivariateKernel h samples x y
--                          / ( univariateKernel h (V.map fst samples) x
--                            * univariateKernel h (V.map snd samples) y
--                            )

-- -- | The joint probability distribution estimation with a univariate standard
-- -- normal distribution.
-- univariateKernel :: Bandwidth -> V.Vector Double -> Double -> Double
-- univariateKernel (Bandwidth h) samples x =
--     (1 / (fromIntegral . V.length $ samples))
--         * (1 / ((sqrt (2 * pi)) * h))
--         * (V.sum . V.map gauss $ samples)
--   where
--     gauss sx     = exp (((\a -> traceShow (x, sx) a) $ (x - sx) ^ 2) / (2 * (h ^ 2)))

-- -- | The joint probability distribution estimation with a bivariate standard
-- -- normal distribution.
-- bivariateKernel :: Bandwidth
--                 -> V.Vector (Double, Double)
--                 -> Double
--                 -> Double
--                 -> Double
-- bivariateKernel (Bandwidth h) samples x y =
--     (1 / (fromIntegral . V.length $ samples))
--         * (1 / (2 * pi * (h ^ 2)))
--         * ( V.sum . V.map (uncurry gauss) $ samples)
--   where
--     gauss sx sy = exp ((((x - sx) ^ 2) + ((y - sy) ^2)) / (2 * (h ^ 2)))

-- -- | The bivariate standard normal density.
-- standardUnivariateGauss :: Double -> V.Vector Double -> Double
-- standardUnivariateGauss x samples = (1 / (sig * sqrt (2 * pi)))
--                                   * exp (- (((x - mu) ^ 2) / (2 * (sig ^ 2))))
--   where
--     mu  = 0 -- mean samples
--     sig = 1 -- stdDev samples

-- -- | The bivariate standard normal density.
-- standardBivariateGauss :: Double
--                        -> Double
--                        -> V.Vector (Double, Double)
--                        -> Double
-- standardBivariateGauss x y samples =
--     (1 / (2 * pi * sig1 * sig2 * (sqrt (1 - (rho ^ 2)))))
--         * exp (-(z / (2 * (1 - (rho ^ 2)))))
--   where
--     z    = (((x - mu1) ^ 2) / (sig1 ^ 2))
--          - ((2 * rho * (x - mu1) * (y - mu2)) / (sig1 * sig2))
--          + (((y - mu2) ^ 2) / (sig2 ^ 2))
--     mu1  = 0 -- mean . fmap fst $ samples
--     mu2  = 0 -- mean . fmap snd $ samples
--     sig1 = 1 -- stdDev . fmap fst $ samples
--     sig2 = 1 -- stdDev . fmap snd $ samples
--     rho  = 0 -- covariance samples
--          -- / ((stdDev . fmap fst $ samples) * (stdDev . fmap snd $ samples))

-- | Complete R version.
-- Take one level and get the similarity matrix by using mutual information (a
-- co-expression network).
-- Does *not* correct for positions with entityDiff here.
getSimMatAracneR :: Size -> StandardLevel -> R.R s EdgeSimMatrix
getSimMatAracneR size level = do
    rDF <- standardLevelToRJSON level

    [r| suppressPackageStartupMessages(library("minet")) |]
    rMat <- [r| df = rDF_hs;
                write("Starting minet.", stderr());
                df = minet(df, method = "aracne", estimator = "mi.shrink", disc = "equalfreq");
                write("Finished minet.", stderr());
                df[is.na(df)] = 0;
                write("Zeros set for minet output.", stderr());
                df
            |]

    res <- rToMat size rMat

    return . EdgeSimMatrix $ res
