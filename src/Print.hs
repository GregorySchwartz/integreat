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
import Control.Monad
import Data.List
import Data.Maybe
import qualified Data.Map.Strict as Map
import Data.Monoid

-- Cabal
import Statistics.Resampling.Bootstrap
import qualified Data.Text as T
import TextShow
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS

-- Local
import Types
import Utility

-- | Print the node correlation scores after integration.
printNodeCorrScores :: IDVec
                    -> Maybe UnifiedData
                    -> NodeCorrScoresMap
                    -> V.Vector NodeCorrScoresInfo
                    -> T.Text
printNodeCorrScores (IDVec idVec) unified scoreMap infos =
    T.append (header <> "\n")
        . T.unlines
        . zipWith
            (\x y -> T.intercalate "," $ x : fmap (fromMaybe "" . fmap showt) y)
            vertexCols
        . fmap getRow
        . V.toList
        $ infos
  where
    header = T.intercalate "," $ startHeader <> endHeader
    startHeader = ("vertex,numSamples" :)
                $ ( fmap ((\(LevelName !x, LevelName !y) -> x <> "_" <> y) . fst)
                  . Map.toAscList
                  . unNodeCorrScoresMap
                  $ scoreMap
                  )
    endHeader   = [ "average"
                  , "estPoint"
                  , "CILower"
                  , "CIUpper"
                  , "CIWidth"
                  , "rankProd"
                  , "pValRankProd"
                  ]
    vertexCols  = fmap vertexCol [0..(V.length idVec - 1)]
    vertexCol x = T.intercalate "," [ unID . (V.!) idVec $ x
                                    , maybe "0" (T.pack . show . Map.size)
                                    . (=<<) ( Map.lookup ((V.!) idVec x)
                                            . unUnifiedData
                                            )
                                    $ unified
                                    ]
    getRow info = (fmap Just . nodeCorrScore $ info)
               <> [ avgNodeCorrScores info
                  , fmap (estPoint . unBootstrap)
                  . avgStatisticNodeCorrScores
                  $ info
                  , ciLower info
                  , ciUpper info
                  , ciWidth info
                  , rankProdNodeCorrScores info
                  , rankProdPValNodeCorrScores info
                  ]
    ciLower info =
        fmap (estLowerBound . unBootstrap) . avgStatisticNodeCorrScores $ info
    ciUpper info =
        fmap (estUpperBound . unBootstrap) . avgStatisticNodeCorrScores $ info
    ciWidth info = fmap abs $ (-) <$> ciLower info <*> ciUpper info
