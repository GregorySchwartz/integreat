{- Integrate
Gregory W. Schwartz

Collections the functions pertaining to the integration of levels of data.
-}

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE QuasiQuotes #-}

module Integrate
    ( integrateCosineSim
    , integrateWalker
    -- , integrateCSRW
    , getNodeCorrScoresInfo
    ) where

-- Standard
import Data.Maybe
import qualified Data.Set as Set
import qualified Data.Map.Strict as Map

-- Cabal
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.Text as T
import qualified Data.Aeson as A
import Data.Graph.Inductive

import qualified Foreign.R as R
import Language.R.Instance as R
import Language.R.Literal as R
import Language.R.QQ

-- Local
import Types
import Utility
import Alignment.Cosine
import Alignment.RandomWalk
-- import Alignment.CSRW

-- | Do integration using cosine similarity.
integrateCosineSim :: Threads
                   -> Permutations
                   -> Size
                   -> VertexSimMap
                   -> EdgeSimMap
                   -> IO NodeCorrScoresMap
integrateCosineSim threads nPerm size vertexSimMap (EdgeSimMap edgeSimMap) =
    fmap (NodeCorrScoresMap . Map.fromList) . parPairsM threads compare $ levels
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
        return ( (x, y), res)
    levels        = Map.keys edgeSimMap
    lookupLevel l = lookupWithError
                        ( "Level: "
                       ++ (T.unpack $ unLevelName l)
                       ++ " notFound in integration step."
                        )
                        l
                        edgeSimMap

-- | Do integration using a random walker.
integrateWalker :: WalkerRestart -> Counter -> GrMap -> IO NodeCorrScoresMap
integrateWalker walkerRestart counterStop (GrMap grMap) =
    fmap (NodeCorrScoresMap . Map.fromList)
        . pairsM eval
        . Map.toAscList
        $ grMap
  where
    eval (!x, gr1) (!y, gr2) =
        fmap ((x, y),)
            . fmap (NodeCorrScores . fmap (, Nothing) . V.fromList)
            . sequence
            . fmap (randomWalkerScore walkerRestart counterStop gr1 gr2)
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

-- | Get the rank product and their p values of node correspondence scores. Uses
-- the RankProd library from R.
getRankProdNodeCorrScores :: [NodeCorrScores]
                          -> R s (FlatNodeCorrScores, PValNodeCorrScores)
getRankProdNodeCorrScores xs = do
    let crate = B.unpack . A.encode . fmap (fmap fst . unNodeCorrScores) $ xs

    resultsDF <- [r| suppressPackageStartupMessages(library("jsonlite"))
                     suppressPackageStartupMessages(library("RankProd"))
                     write("Getting rank product.", stderr());
                     df = fromJSON(crate_hs)
                     x = capture.output(result <- RP(t(df), rep(1, ncol(t(df))), logged=TRUE))
                     result
                 |]

    flat  <- [r| resultsDF_hs[["RPs"]][,1] |]
    pVals <- [r| resultsDF_hs[["pval"]][,1] |]

    return
        ( FlatNodeCorrScores . V.fromList $ (R.fromSomeSEXP flat :: [Double])
        , PValNodeCorrScores . V.fromList $ (R.fromSomeSEXP pVals :: [Double])
        )

-- | Get all information from the node correspondence scores.
getNodeCorrScoresInfo :: NodeCorrScoresMap -> R s NodeCorrScoresInfo
getNodeCorrScoresInfo scoreMap = do
    let avgScores =
            getAvgNodeCorrScores
                . Map.elems
                . unNodeCorrScoresMap
                $ scoreMap
        avgPVals =
            PValNodeCorrScores
                . avgVecVec
                . fmap ( fmap (fromMaybe 1 . fmap unPValue . snd)
                       . unNodeCorrScores
                       )
                . Map.elems
                . unNodeCorrScoresMap
                $ scoreMap

    (rankProds,pVals) <-
        getRankProdNodeCorrScores . Map.elems . unNodeCorrScoresMap $ scoreMap

    return $
        NodeCorrScoresInfo
        { nodeCorrScoresMap = scoreMap
        , avgNodeCorrScores = avgScores
        , avgPValNodeCorrScores = avgPVals
        , rankProdNodeCorrScores = rankProds
        , rankProdPValNodeCorrScores = pVals
        }
