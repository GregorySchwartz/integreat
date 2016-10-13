{- Print
Gregory W. Schwartz

Collections the functions pertaining to the printing of results.
-}

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE OverloadedStrings #-}

module Print
    ( printNodeCorrScores
    ) where

-- Standard
import Data.List
import qualified Data.Map.Strict as Map
import Data.Monoid

-- Cabal
import qualified Data.Text as T
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import TextShow

-- Local
import Types
import Utility

-- | Print the node correlation scores after integration.
printNodeCorrScores :: IDVec -> NodeCorrScoresInfo -> T.Text
printNodeCorrScores (IDVec idVec) info =
    T.append (header <> "\n")
        . T.unlines
        . foldl'
            (\acc -> zipWith (\x y -> x <> "," <> showt y) acc . V.toList)
            (fmap (unID . (V.!) idVec) [0..(V.length idVec - 1)])
        $ cols
  where
    header = T.intercalate ","
           . ("vertex" :)
           $ ( fmap ((\(LevelName !x, LevelName !y) -> x <> "_" <> y) . fst)
             . Map.toAscList
             . unNodeCorrScoresMap
             . nodeCorrScoresMap
             $ info
             )
          <> ["average", "rankProd", "pValRankProd"]
    cols = ( fmap (unNodeCorrScores . snd)
           . Map.toAscList
           . unNodeCorrScoresMap
           . nodeCorrScoresMap
           $ info
           )
        <> [ unFlatNodeCorrScores . avgNodeCorrScores $ info
           , unFlatNodeCorrScores . rankProdNodeCorrScores $ info
           , unPValNodeCorrScores . rankProdPValNodeCorrScores $ info
           ]
