{- Alignment.RandomWalk
Gregory W. Schwartz

Collections the functions pertaining to random walk over the network.
-}

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DuplicateRecordFields #-}

module Alignment.RandomWalk
    ( randomWalkerScore
    ) where

-- Standard
import Data.Maybe
import qualified Data.Set as Set
import Data.List
import Data.Tuple
import qualified Data.IntMap.Strict as IMap

-- Cabal
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Data.Text as T
import Data.Graph.Inductive
import Data.Random
import Control.Monad.Reader
import Control.Monad.State
import Control.Lens

-- Local
import Types
import Utility

-- | Update the current position of the walker and the distribution.
updateWalker :: Int -> Walker ()
updateWalker v = do
    oldState <- get

    let addOne Nothing  = Just 0
        addOne (Just x) = Just $ x + 1

    put
        . WalkerState
        . over _2 (IMap.alter addOne v)
        . set _1 v
        . unWalkerState
        $ oldState

-- | Transform edge weights to probability thresholds, setting the last element
-- to 0.
weightsToProbs :: Adj Double -> Adj Double
weightsToProbs xs = flip zip vertices
                  . setLastZero
                  . drop 1
                  . scanl' (\acc x -> acc - x) 1
                  . fmap (/ sum normalizedWeights)
                  $ normalizedWeights
  where
    setLastZero = reverse . (0 :) . drop 1 . reverse
    normalizedWeights = minMaxNorm . fmap fst $ xs
    vertices    = fmap snd xs

-- | Get a random neighbor based on their edge weights.
randomNeighbor :: Adj Double -> IO (Maybe Int)
randomNeighbor [] = return Nothing
randomNeighbor ns = do
    stepQuery <- sample (Uniform 0 1 :: Uniform Double)
    return
        . Just
        . snd
        . head
        . dropWhile ((>= stepQuery) . fst)
        . weightsToProbs
        $ ns

-- | Get the distribution of where the walker visits.
getWalkDist :: Counter -> Counter -> Walker ()
getWalkDist counterStop counter = do
    gr <- fmap (unLevelGr . eGr) ask
    v  <- fmap (view _1 . unWalkerState) get

    restartThreshold <-
        fmap (unWalkerRestart . (restart :: (Environment -> WalkerRestart))) ask
    restartQuery     <- liftIO $ sample (Uniform 0 1 :: Uniform Double)

    v' <- if restartQuery <= restartThreshold
            then do
                fmap v0 ask
            else do
                liftIO . fmap (fromMaybe v) . randomNeighbor . lneighbors gr $ v

    updateWalker v'

    unless (counter > counterStop)
        . getWalkDist counterStop
        $ (counter + Counter 1)

-- | Run the random walker on a vertex.
evalWalker :: WalkerRestart
           -> Counter
           -> LevelGr
           -> Int
           -> IO (IMap.IntMap Int)
evalWalker walkerRestart counterStop gr v = do
    distribution <- fmap (view _2 . unWalkerState)
                  . flip execStateT (st v)
                  . runReaderT
                        (unWalker $ getWalkDist counterStop (Counter 0))
                  $ env
    return distribution
  where
    st v = WalkerState (v, IMap.empty)
    env  = Environment { eGr     = gr
                       , restart = walkerRestart
                       , v0      = v
                       }

-- | Get the node correspondence score of a vertex between two levels.
randomWalkerScore :: WalkerRestart
                  -> Counter
                  -> LevelGr
                  -> LevelGr
                  -> Int
                  -> IO Double
randomWalkerScore walkerRestart counterStop gr1 gr2 v = do
    gr1Dist <- fmap (IMap.map fromIntegral)
             . evalWalker walkerRestart counterStop gr1
             $ v
    gr2Dist <- fmap (IMap.map fromIntegral)
             . evalWalker walkerRestart counterStop gr2
             $ v

    return . cosineSimIMap gr1Dist $ gr2Dist
