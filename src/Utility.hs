{- Utility
Gregory W. Schwartz

Collections all miscellaneous functions.
-}

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes #-}

module Utility
    ( minMaxNorm
    , subsetIMap
    , imapToVec
    , imapMatToMat
    , lookupWithError
    , largestLeftEig
    , (/.)
    , removeMatchFilter
    , applyRows
    , avgVec
    , avgVecVec
    , getVertexSim
    , pairs
    , pairsM
    , parPairsM
    , triples
    , flipToo
    , listToTuple
    , sameWithEntityDiff
    , groupDataSets
    , fillIntMap
    , getAccuracy
    , nub'
    ) where

-- Standard
import Data.Maybe
import Data.List
import qualified Data.Set as Set
import qualified Data.Map.Strict as Map
import qualified Data.IntMap.Strict as IMap
import Data.Function (on)
import Control.Concurrent.Async

-- Cabal
import qualified Data.Vector as V
import qualified Data.Vector.Generic as VG
import qualified Data.Vector.Storable as VS
import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.Text as T
import Data.ByteString.Lazy.Search (replace)
import Data.Graph.Inductive
import qualified Data.Aeson as JSON
import Control.Lens
import Numeric.LinearAlgebra

-- Local
import Types

-- | Min max normalize.
minMaxNorm :: [Double] -> [Double]
minMaxNorm xs = fmap (\x -> (x - minimum xs) / (maximum xs - minimum xs)) xs

-- | Get a subset of indices of an IntMap.
subsetIMap :: [Int] -> IMap.IntMap a -> IMap.IntMap a
subsetIMap idx = flip IMap.intersection (IMap.fromList . zip idx $ [0..])

-- | Convert an IntMap to a vector with a filler.
imapToVec :: (VG.Vector v a, Show a) => Size -> a -> IMap.IntMap a -> v a
imapToVec (Size size) filler =
    (VG.//) (VG.replicate size filler) . IMap.toAscList
    
-- | Convert an IntMap of IntMap to a square matrix with a filler.
imapMatToMat :: Size
             -> Double
             -> IMap.IntMap (IMap.IntMap Double)
             -> Matrix Double
imapMatToMat (Size size) filler =
    assoc (size, size) filler
        . concatMap (\ (!row, !v)
                    -> fmap (\ (!col, !v)
                            -> ((row, col), v)
                            ) . IMap.toAscList $ v
                    )
        . IMap.toAscList

-- | Map lookup with a custom error if the value is not found.
lookupWithError :: (Ord a) => String -> a -> Map.Map a b -> b
lookupWithError err x = fromMaybe (error err) . Map.lookup x

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

-- | Remove indices matching a boolean function from both vectors but make
-- sure that the indices match. NOT NEEDED WITH COSINE.
removeMatchFilter :: (VS.Storable a)
                  => (a -> Bool)
                  -> Vector a
                  -> Vector a
                  -> (Vector a, Vector a)
removeMatchFilter f xs = over _2 VS.convert
                       . over _1 VS.convert
                       . V.unzip
                       . filterBad (V.convert xs)
                       . V.convert
  where
    filterBad x = V.filter (\(!a, !b) -> (not . f $ a) && (not . f $ b))
                . V.zip x

-- | Apply a folding function to a list of row vectors.
applyRows :: (Element a, VS.Storable b)
          => (Vector a -> b)
          -> [Vector a]
          -> Vector b
applyRows f = fromList . fmap f . toColumns . fromRows

-- | Average of a vector.
avgVec :: (Num a, Fractional a, VG.Vector v a) => v a -> a
avgVec xs = VG.sum xs / (fromIntegral $ VG.length xs)

-- | Average entries of a list of vectors.
avgVecVec :: (Num a, Fractional a, VG.Vector v a) => [v a] -> v a
avgVecVec xs = VG.map (/ genericLength xs)
             . foldl1' (VG.zipWith (+))
             $ xs

-- | Get the vertex similarity matrix for two levels, erroring out if the
-- levels don't exist.
getVertexSim :: LevelName -> LevelName -> VertexSimMap -> VertexSimValues
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

-- | From
-- http://stackoverflow.com/questions/34044366/how-to-extract-all-unique-pairs-of-a-list-in-haskell,
-- extract the unique pairings of a list and apply a function to them.
pairs :: (a -> a -> b) -> [a] -> [b]
pairs f l = [f x y | (x:ys) <- tails l, y <- ys]

-- | Extract the unique pairings of a list and apply a function to them within a
-- monad.
pairsM :: (Monad m) => (a -> a -> m b) -> [a] -> m [b]
pairsM f l = sequence [f x y | (x:ys) <- tails l, y <- ys]

-- | Extract the unique pairings of a list and apply a function to them within a
-- monad. Parallel version.
parPairsM :: (a -> a -> IO b) -> [a] -> IO [b]
parPairsM _ [] = return []
parPairsM f l  = mapConcurrently (uncurry evalM) . concatMap go . tails $ l
  where
    go [] = []
    go [_] = []
    go (x:ys) = fmap (x,) ys
    evalM !x !y = do
        !res <- f x y
        return res

-- | Extract the unique triplets of a list and apply a function to them.
triples :: (a -> a -> a -> b) -> [a] -> [b]
triples f l = [f x y z | (x:ys) <- tails l, (y:zs) <- tails ys, z <- zs]

-- | Take a tuple index with a value and return it with its flip.
flipToo :: ((a, a), b) -> [((a, a), b)]
flipToo all@((!x, !y), !z) = [all, ((y, x), z)]

-- | Convert a list to a tuple.
listToTuple :: (Show a) => [a] -> (a, a)
listToTuple [!x, !y] = (x, y)
listToTuple (x:_)    = error ("Wrong pairing for " ++ show x)

-- | Check if two entities are actually the same if one contains the entityDiff
-- while the other does not.
sameWithEntityDiff :: Maybe EntityDiff -> ID -> ID -> Bool
sameWithEntityDiff Nothing (ID e1) (ID e2)                           = False
sameWithEntityDiff (Just (EntityDiff eDiff)) (ID e1) (ID e2)
    | T.count eDiff e1 == 0 && T.count eDiff e2 == 0                 = False
    | T.count eDiff e1 > 0 && T.count eDiff e2 > 0                   = False
    | (head . T.splitOn eDiff $ e1) == (head . T.splitOn eDiff $ e2) = True
    | otherwise                                                      = False

-- | Group together data sets.
groupDataSets :: Maybe Entity -> Maybe Entity -> Maybe (Double, Double)
groupDataSets Nothing _           = Nothing
groupDataSets _ Nothing           = Nothing
groupDataSets (Just !x) (Just !y) = Just (_entityValue x, _entityValue y)

-- | Fill in an intmap with zeros for missing values.
fillIntMap :: Size
           -> IMap.IntMap (Double, Maybe Statistic)
           -> IMap.IntMap (Double, Maybe Statistic)
fillIntMap (Size size) m =
    foldl'
        (\acc x -> IMap.alter (maybe (Just (0, Nothing)) Just) x acc)
        m
        [0..size - 1]

-- | Rank the node correspondence scores.
rankNodeCorrScores :: IDVec -> NodeCorrScores -> [(Double, ID)]
rankNodeCorrScores (IDVec idVec) = zip [1..]
                                 . fmap fst
                                 . sortBy (compare `on` snd)
                                 . V.toList
                                 . V.imap (\ !i !v -> (idVec V.! i, v))
                                 . VS.convert
                                 . fmap fst
                                 . unNodeCorrScores

-- | Get the accuracy of a run. In this case, we get the total rank below the
-- number of permuted vertices divided by the theoretical maximum (so if there were
-- five changed vertices out off 10 and two were rank 8 and 10 while the others
-- were in the top five, we would have (1 - ((3 + 5) / (10 + 9 + 8 + 7 + 6))) as
-- the accuracy."
getAccuracy :: Set.Set ID -> IDVec -> V.Vector NodeCorrScoresInfo -> Double
getAccuracy truth (IDVec idVec) =
    (1 -)
        . (/ fact)
        . sum
        . filter (> 0)
        . fmap ( flip
                    (-)
                    (fromIntegral . Set.size $ truth)
               . fst
               )
        . filter (flip Set.member truth . snd)
        . rankNodeCorrScores (IDVec idVec)
        . NodeCorrScores
        . fmap ((, Nothing) . fromJust . avgNodeCorrScores)
  where
    fact = fromIntegral
         . sum
         . fmap ((V.length idVec - length truth) -)
         $ [0..(length truth - 1)]

-- | Faster nub.
nub' :: (Eq a, Ord a) => [a] -> [a]
nub' = Set.toList . Set.fromList
