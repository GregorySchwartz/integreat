{- NetworkGeneration
Gregory W. Schwartz

Collections the functions pertaining to the generation of vertices and
edges in the data sets using a similarity matrix.
-}

{-# LANGUAGE BangPatterns #-}

module NetworkGeneration
    ( getSimMat
    ) where

-- Standard
import qualified Data.Map.Strict as Map
import qualified Data.Sequence as Seq
import qualified Data.Foldable as F

-- Cabal
import qualified Data.Vector as V
import Statistics.Correlation.Kendall
import Numeric.LinearAlgebra

-- Local
import Types
import Utility

-- | Take one level and get the similarity matrix by using correlations (a
-- co-expression network). The default value is applied to missing data.
-- Consider using a value that would never occur: for instance with the
-- Kendall correlation, -5 would never appear.
getSimMat :: Default -> IDMap -> Level -> EdgeSimMatrix
getSimMat (Default def) (IDMap idMap) (Level level) =
    EdgeSimMatrix
        . assoc (Map.size idMap, Map.size idMap) def
        $ getCorr <$> Map.toList level <*> Map.toList level
  where
      getCorr (k1, x) (k2, y) = ( ( lookupWithError (keyNotFound k1) k1 idMap
                                  , lookupWithError (keyNotFound k2) k2 idMap
                                  )
                                  , correlate x y
                                )
      keyNotFound k = "ID: " ++ show k ++ " not found."

-- | Correlate two groups of entities, where each group is a collection of
-- measurements (specifically for a single type of entity, for instance
-- a single protein).
correlate :: Seq.Seq Entity -> Seq.Seq Entity -> Double
correlate e1 e2 =
    kendall
        . V.zip (V.fromList . F.toList . fmap _entityValue $ e1)
        . V.fromList
        . F.toList
        . fmap _entityValue
        $ e2
