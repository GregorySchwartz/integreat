{- Alignment.Cosine
Gregory W. Schwartz

Collections the functions pertaining to getting integration by cosine
similarity.
-}

{-# LANGUAGE BangPatterns #-}

module Alignment.Cosine
    ( cosineIntegrate
    , cosineSim
    , cosineSimIMap
    , cosinePerm
    ) where

-- Standard
import Data.Tuple
import Data.List
import qualified Data.IntMap.Strict as IMap
import Control.Monad

-- Cabal
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import Data.Random
import Control.Lens
import Numeric.LinearAlgebra

-- Local
import Types
import Utility

-- | Get the cosine similarity of two levels.
cosineIntegrate :: Permutations
                -> Size
                -> VertexSimMap
                -> LevelName
                -> LevelName
                -> EdgeSimMatrix
                -> EdgeSimMatrix
                -> IO NodeCorrScores
cosineIntegrate nPerm size vMap l1 l2 e1 e2 =
    fmap NodeCorrScores
        . sequence
        . V.fromList
        . fmap snd
        . IMap.toAscList
        . fillIntMap size
        . IMap.intersectionWith (cosinePerm nPerm edgeVals) newE1
        $ newE2
  where
    edgeVals :: EdgeValues
    edgeVals = EdgeValues . concatMap IMap.elems . IMap.elems $ newE2
    newE1 :: IMap.IntMap (IMap.IntMap Double)
    newE1         =
        unEdgeSimMatrix . cosineUpdateSimMat e1 . unVertexSimValues $ vertexSim
    newE2 :: IMap.IntMap (IMap.IntMap Double)
    newE2         =
        unEdgeSimMatrix . cosineUpdateSimMat e2 . unVertexSimValues $ vertexSim
    vertexSim     = getVertexSim l1 l2 vMap

-- | Update the similarity matrices where the diagonal contains the similarity
-- between vertices of different levels.
cosineUpdateSimMat :: EdgeSimMatrix -> [((Int, Int), Double)] -> EdgeSimMatrix
cosineUpdateSimMat (EdgeSimMatrix edgeSimMat) =
    EdgeSimMatrix
        . foldl' (\acc ((!x, !y), !z) -> IMap.alter (f y z) x acc) edgeSimMat
        . (\xs -> xs ++ fmap (over _1 swap) xs)
  where
    f y z Nothing  = Just . IMap.singleton y $ z
    f y z (Just v) = Just . IMap.insert y z $ v

-- | Cosine similarity.
cosineSim :: Vector Double -> Vector Double -> Double
cosineSim x y = dot x y / (norm_2 x * norm_2 y)

-- | Cosine similarity of two IntMaps.
cosineSimIMap :: IMap.IntMap Double -> IMap.IntMap Double -> Double
cosineSimIMap x y = (imapSum $ IMap.intersectionWith (*) x y)
                  / (imapNorm x * imapNorm y)
  where
    imapNorm = sqrt . imapSum . IMap.map (^ 2)
    imapSum  = IMap.foldl' (+) 0

-- | Get the cosine similarity and the p value of that similarity through the
-- permutation test.
cosinePerm :: Permutations
           -> EdgeValues
           -> IMap.IntMap Double
           -> IMap.IntMap Double
           -> IO (Double, Maybe PValue)
cosinePerm (Permutations nPerm) edgeVals xs ys = do
    let obs     = cosineSimIMap xs ys
        expTest = (>= abs obs) . abs

    -- Lotsa space version.
    -- vals <- mapM (const (shuffleCosine xs ys)) . V.replicate nPerm $ 0
    -- let exps = V.filter (\x -> abs x >= abs obs) vals

    let successes :: Int -> Int -> IO Int
        successes !acc 0  = return acc
        successes !acc !n = do
            res <- randomCosine edgeVals xs $ ys
            if expTest res
                then successes (acc + 1) (n - 1)
                else successes acc (n - 1)

    exp <- successes 0 nPerm

    let pVal = PValue $ (fromIntegral exp) / (fromIntegral nPerm)

    return (obs, Just pVal)

-- | Random sample cosine.
randomCosine :: EdgeValues
             -> IMap.IntMap Double
             -> IMap.IntMap Double
             -> IO Double
randomCosine (EdgeValues edgeVals) xs ys = do
    shuffledYS <- getSampleMap edgeVals ys

    return . cosineSimIMap xs $ shuffledYS

-- | Sample map values from the complete set of values.
getSampleMap :: [a] -> IMap.IntMap a -> IO (IMap.IntMap a)
getSampleMap vals m = do
    let keys = IMap.keys m

    newElems <- sample . replicateM (IMap.size m) . randomElement $ vals

    return . IMap.fromList . zip keys $ newElems
    
-- | Random shuffle cosine.
shuffleCosine :: IMap.IntMap Double
              -> IMap.IntMap Double
              -> IO Double
shuffleCosine xs ys = do
    shuffledXS <- shuffleMap xs
    shuffledYS <- shuffleMap ys

    return . cosineSimIMap shuffledXS $ shuffledYS

-- | Shuffle map keys.
shuffleMap :: IMap.IntMap a -> IO (IMap.IntMap a)
shuffleMap m = do
    let (keys, elems) = unzip . IMap.toList $ m

    newKeys <- sample . shuffle $ keys

    return . IMap.fromList . zip newKeys $ elems
