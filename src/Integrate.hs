{- Integrate
Gregory W. Schwartz

Collections the functions pertaining to the integration of levels of data.
-}

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE PackageImports #-}

module Integrate
    ( integrateCosineSim
    , integrateWalker
    , integrateWalkerSim
    -- , integrateCSRW
    , getNodeCorrScoresInfo
    ) where

-- Standard
import Data.List
import Data.Maybe
import Data.Monoid
import qualified Data.Map.Strict as Map
import qualified Data.Set as Set

-- Cabal
import Data.Graph.Inductive
import qualified Statistics.RankProduct as RP
import qualified "rank-product" Statistics.Types as RP
import qualified Data.Aeson as A
import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.Text as T
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS

-- import qualified Foreign.R as R
-- import Language.R.Instance as R
-- import Language.R.Literal as R
-- import Language.R.QQ

-- Local
import Types
import Utility
import Alignment.Cosine
import Alignment.RandomWalk
-- import Alignment.CSRW

-- | Do integration using cosine similarity.
integrateCosineSim :: Permutations
                   -> Size
                   -> VertexSimMap
                   -> EdgeSimMap
                   -> IO NodeCorrScoresMap
integrateCosineSim nPerm size vertexSimMap (EdgeSimMap edgeSimMap) =
    fmap (NodeCorrScoresMap . Map.fromList) . pairsM compare $ levels
  where
    compare !x !y = do
        res <- cosineIntegrate
                nPerm
                size
                vertexSimMap
                x
                y
                (lookupLevel x)
                (lookupLevel y)
        return ((x, y), res)
    levels        = Map.keys edgeSimMap
    lookupLevel l = lookupWithError
                        ( "Level: "
                       ++ (T.unpack $ unLevelName l)
                       ++ " notFound in integration step."
                        )
                        l
                        edgeSimMap

-- | Do integration using a random walker, HotNet2 style.
integrateWalker :: Permutations
                -> Size
                -> WalkerRestart
                -> EdgeSimMap
                -> IO NodeCorrScoresMap
integrateWalker nPerm size walkerRestart (EdgeSimMap edgeSimMap) =
    fmap (NodeCorrScoresMap . Map.fromList) . pairsM compare $ levels
  where
    compare !x !y = do
        res <- randomWalkerScore
                nPerm
                size
                walkerRestart
                (lookupLevel x)
                (lookupLevel y)
        return ((x, y), res)
    levels        = Map.keys edgeSimMap
    lookupLevel l = lookupWithError
                        ( "Level: "
                       ++ (T.unpack $ unLevelName l)
                       ++ " notFound in integration step."
                        )
                        l
                        edgeSimMap

-- | Do integration using a random walker, simulated.
integrateWalkerSim :: WalkerRestart -> Counter -> GrMap -> IO NodeCorrScoresMap
integrateWalkerSim walkerRestart counterStop (GrMap grMap) =
    fmap (NodeCorrScoresMap . Map.fromList)
        . pairsM eval
        . Map.toAscList
        $ grMap
  where
    eval (!x, gr1) (!y, gr2) =
        fmap ((x, y),)
            . fmap (NodeCorrScores . fmap (, Nothing) . V.fromList)
            . sequence
            . fmap (randomWalkerScoreSim walkerRestart counterStop gr1 gr2)
            . allNodes (unLevelGr gr1)
            . unLevelGr
            $ gr2
    allNodes gr1 gr2 = Set.toList
                     . Set.union (Set.fromList . nodes $ gr2)
                     . Set.fromList
                     . nodes
                     $ gr1

-- -- | Use integration using a random walker "CSRW".
-- integrateCSRW
--     :: VertexSimMap
--     -> EdgeSimMap
--     -> WalkerRestart
--     -> Counter
--     -> IO NodeCorrScoresMap
-- integrateCSRW vertexSimMap edgeSimMap restart counterStop =
--       fmap (NodeCorrScoresMap . Map.fromList) . pairsM eval $ levels
--   where
--     eval !l1 !l2   = fmap ((l1, l2),)
--                    . evalWalker vertexSimMap edgeSimMap restart counterStop l1
--                    $ l2
--     levels         = Map.keys . unEdgeSimMap $ edgeSimMap
--     lookupLevel l  = lookupWithError
--                         ("Level: " ++ (T.unpack $ unLevelName l) ++ " notFound")
--                         l
--                    . unEdgeSimMap
--                    $ edgeSimMap

-- | Get the average node correspondence scores.
getAvgNodeCorrScores :: [NodeCorrScores] -> FlatNodeCorrScores
getAvgNodeCorrScores =
    FlatNodeCorrScores . avgVecVec . fmap (fmap fst . unNodeCorrScores)

-- | Get the rank product and their p values of node correspondence scores.
getRankProdNodeCorrScores :: [NodeCorrScores]
                          -> IO (FlatNodeCorrScores, PValNodeCorrScores)
getRankProdNodeCorrScores xs = do
    res <- RP.rankProductPermutation (RP.Permutations 1000)
         . fmap RP.Entity
         . transpose
         . fmap (V.toList . fmap fst . unNodeCorrScores)
         $ xs

    let flat  = fmap (RP.unRankProductEntity . fst) res
        pVals = fmap (RP.unPValue . snd) res

    return
        ( FlatNodeCorrScores . V.fromList $ flat
        , PValNodeCorrScores . V.fromList $ pVals
        )

-- -- | Get the rank product and their p values of node correspondence scores. Uses
-- -- the RankProd library from R.
-- getRankProdNodeCorrScoresR :: [NodeCorrScores]
--                            -> R s (FlatNodeCorrScores, PValNodeCorrScores)
-- getRankProdNodeCorrScoresR xs = do
--     let crate = B.unpack . A.encode . fmap (fmap fst . unNodeCorrScores) $ xs

--     resultsDF <- [r| suppressPackageStartupMessages(library("jsonlite"));
--                      suppressPackageStartupMessages(library("RankProd"));
--                      suppressPackageStartupMessages(library("data.table"));
--                      write("Getting rank product.", stderr());
--                      crate = copy(crate_hs);
--                      df = fromJSON(crate);
--                      x = capture.output(result <- RP(t(df), rep(1, ncol(t(df))), logged=TRUE));
--                      result
--                  |]

--     flat  <- [r| resultsDF_hs[["RPs"]][,1] |]
--     pVals <- [r| resultsDF_hs[["pval"]][,1] |]

--     return
--         ( FlatNodeCorrScores . V.fromList $ (R.fromSomeSEXP flat :: [Double])
--         , PValNodeCorrScores . V.fromList $ (R.fromSomeSEXP pVals :: [Double])
--         )

-- | Get all information from the node correspondence scores.
getNodeCorrScoresInfo :: NodeCorrScoresMap
                      -> IO (V.Vector NodeCorrScoresInfo)
getNodeCorrScoresInfo scoreMap = do
    let scores = Map.toAscList . unNodeCorrScoresMap $ scoreMap
        avgScores = getAvgNodeCorrScores . fmap snd $ scores
        -- avgStatistic =
        --     StatisticNodeCorrScores
        --         . avgVecVec
        --         . fmap ( fmap (fromMaybe 0 . fmap unPValue . snd)
        --                . unNodeCorrScores
        --                )
        --         . fmap snd
        --         $ scores
        avgStatistic :: StatisticNodeCorrScores
        avgStatistic =
            StatisticNodeCorrScores
                . fmap snd
                . unNodeCorrScores
                . snd
                . head
                $ scores
        nodeScores = fmap (fmap fst . unNodeCorrScores . snd) scores

    (rankProds, pVals) <- getRankProdNodeCorrScores . fmap snd $ scores

    let f avgScoresX avgStatisticX rankProdsX pValsX rest =
            NodeCorrScoresInfo
                { nodeCorrScore              = rest
                , avgNodeCorrScores          = Just avgScoresX
                , avgStatisticNodeCorrScores = avgStatisticX
                , rankProdNodeCorrScores     = Just rankProdsX
                , rankProdPValNodeCorrScores = Just pValsX
                }
        nodeCorrScoresInfo =
            V.zipWith5 f
                    (unFlatNodeCorrScores avgScores)
                    (unStatisticNodeCorrScores avgStatistic)
                    (unFlatNodeCorrScores rankProds)
                    (unPValNodeCorrScores pVals)
                . V.fromList
                . transpose
                . fmap V.toList
                $ nodeScores

    return nodeCorrScoresInfo
