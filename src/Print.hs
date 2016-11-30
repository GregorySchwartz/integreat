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
import Control.Monad

-- Cabal
import qualified Data.Text as T
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import TextShow

-- Local
import Types
import Utility

-- | Print the node correlation scores after integration.
printNodeCorrScores :: IDVec
                    -> Maybe UnifiedData
                    -> NodeCorrScoresInfo
                    -> T.Text
printNodeCorrScores (IDVec idVec) unified info =
    T.append (header <> "\n")
        . T.unlines
        . foldl'
            (\ acc -> zipWith (\x y -> x <> "," <> showt y) acc
                    . V.toList
            )
            vertexCols
        $ cols
  where
    header = T.intercalate ","
           . ("vertex,numSamples" :)
           $ ( fmap ((\(LevelName !x, LevelName !y) -> x <> "_" <> y) . fst)
             . Map.toAscList
             . unNodeCorrScoresMap
             . nodeCorrScoresMap
             $ info
             )
          <> ["average", "pValAvg", "rankProd", "pValRankProd"]
    vertexCols  = fmap vertexCol [0..(V.length idVec - 1)]
    vertexCol x = T.intercalate "," [ unID . (V.!) idVec $ x
                                    , maybe "0" (T.pack . show . Map.size)
                                    . (=<<) ( Map.lookup ((V.!) idVec x)
                                            . unUnifiedData
                                            )
                                    $ unified
                                    ]
    cols        = ( fmap (fmap fst . unNodeCorrScores . snd)
                  . Map.toAscList
                  . unNodeCorrScoresMap
                  . nodeCorrScoresMap
                  $ info
                  )
               <> [ unFlatNodeCorrScores . avgNodeCorrScores $ info
                  , unPValNodeCorrScores . avgPValNodeCorrScores $ info
                  , unFlatNodeCorrScores . rankProdNodeCorrScores $ info
                  , unPValNodeCorrScores . rankProdPValNodeCorrScores $ info
                  ]
