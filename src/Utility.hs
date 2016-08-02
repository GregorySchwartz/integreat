{- Utility
Gregory W. Schwartz

Collections all miscellaneous functions.
-}

{-# LANGUAGE BangPatterns #-}

module Utility
    ( lookupWithError
    , getNeighbors
    , largestLeftEig
    , (/.)
    , cosineSim
    , removeMatchFilter
    , applyRows
    , avgVec
    , getVertexSim
    , pairs
    , flipToo
    ) where

-- Standard
import Data.Maybe
import Data.List
import qualified Data.Set as Set
import qualified Data.Map.Strict as Map
import qualified Data.Sequence as Seq

-- Cabal
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Data.Text as T
import Numeric.LinearAlgebra

-- Local
import Types

-- | Map lookup with a custom error if the value is not found.
lookupWithError :: (Ord a) => String -> a -> Map.Map a b -> b
lookupWithError err x = fromMaybe (error err) . Map.lookup x

-- | Get the neighbors of a vertex in the SimilarityMatrix.
getNeighbors :: Int -> EdgeSimMatrix -> Set.Set Int
getNeighbors idx = Set.fromList
                 . V.toList
                 . V.map fst
                 . V.filter ((> 0) . snd)
                 . V.imap (,)
                 . VS.convert
                 . flip (!) idx
                 . unEdgeSimMatrix

-- | Get the largest left eigenvector from an eig funciton call. The matrix, a
-- transition probability matrix in this program, is assumed to be symmetrical
-- here.
largestLeftEig :: (VS.Vector (Complex Double), Matrix (Complex Double))
           -> VS.Vector Double
largestLeftEig (!eigVal, !eigVec) =
    fst . fromComplex $ (tr eigVec ! maxIndex eigVal)

-- | A more generic division.
(/.) :: (Real a, Fractional c) => a -> a -> c
(/.) x y = fromRational $ toRational x / toRational y

-- | Cosine similarity.
cosineSim :: Vector Double -> Vector Double -> Double
cosineSim x y = dot x y / (norm_2 x * norm_2 y)

-- | Remove indices matching a boolean function from both vectors but make
-- sure that the indices match.
removeMatchFilter :: (VS.Storable a)
                  => (a -> Bool)
                  -> Vector a
                  -> Vector a
                  -> (Vector a, Vector a)
removeMatchFilter !f !x !y = ( filterIdx x
                             , filterIdx y
                             )
  where
    filterIdx = VS.ifilter (\i v -> not . Set.member i $ badIdx)
    badIdx    = Set.fromList
              . VS.toList
              . mappend (VS.findIndices f x)
              . VS.findIndices f $ y

-- | Apply a folding function to a list of row vectors.
applyRows :: (Element a, VS.Storable b)
          => (Vector a -> b)
          -> [Vector a]
          -> Vector b
applyRows f = fromList . fmap f . toColumns . fromRows

-- | Average of a vector.
avgVec :: Vector Double -> Double
avgVec xs = VS.sum xs / (fromIntegral $ VS.length xs)

-- | Get the vertex similarity matrix for two levels, erroring out if the
-- levels don't exist.
getVertexSim :: LevelName -> LevelName -> VertexSimMap -> VertexSimMatrix
getVertexSim l1 l2 (VertexSimMap vMap) =
    case Map.lookup l1 vMap of
        Nothing  -> fromMaybe (error $ levelErr l1)
                  . maybe (error $ levelErr l2) (Map.lookup l1)
                  . Map.lookup l2
                  $ vMap
        (Just x) -> lookupWithError (error $ levelErr l2) l2
                  $ x
  where
    levelErr l    = ( "Level: "
                   ++ (T.unpack $ unLevelName l)
                   ++ " not found in creating vertex similarity map."
                    )

-- | From http://stackoverflow.com/questions/34044366/how-to-extract-all-unique-pairs-of-a-list-in-haskell, extract the unique pairings of a list and apply a function to them.
pairs :: (a -> a -> b) -> [a] -> [b]
pairs f l = [f x y | (x:ys) <- tails l, y <- ys]

-- | Take a tuple index with a value and return it with its flip.
flipToo :: ((a, a), b) -> [((a, a), b)]
flipToo all@((!x, !y), !z) = [all, ((y, x), z)]
