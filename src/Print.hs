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

-- Cabal
import qualified Data.Text as T
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import TextShow

-- Local
import Types
import Utility

-- | Print the node correlation scores after integration.
printNodeCorrScores :: IDVec -> NodeCorrScores -> T.Text
printNodeCorrScores (IDVec idVec) =
    T.append "vertex,value\n"
        . T.unlines
        . V.toList
        . V.imap (\i v -> T.concat [unID . (V.!) idVec $ i, ",", showt v])
        . VS.convert
        . unNodeCorrScores
