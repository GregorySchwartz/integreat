{- Alignment.CSRW
Gregory W. Schwartz

Collections the functions pertaining to random walk over the network.
-}

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DuplicateRecordFields #-}

module Alignment.CSRW
    ( getTransProbMat
    , evalWalker
    , getWalkerNodeCorrespondenceScores
    ) where

-- Standard
import qualified Data.Set as Set
import Data.Tuple

-- Cabal
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Data.Text as T
import qualified Data.Random as R
import Control.Monad.Reader
import Control.Monad.State
import Control.Lens
import Numeric.LinearAlgebra

-- Local
import Types
import Utility

-- | Get the transition probability.
getTransProb :: WalkerChoice -> Set.Set Int -> WalkerCSRW Double
getTransProb Same neighbors = do
    c1   <- fmap (view _1 . unWalkerStateCSRW) get
    vMat <- fmap (unVertexSimMatrix . vertexMat) ask

    return $ vMat ! c1 ! c1 /. (fromIntegral $ Set.size neighbors)
getTransProb choice neighbors = do
    m1 <- fmap (unEdgeSimMatrix . edgeMat1) ask
    m2 <- fmap (unEdgeSimMatrix . edgeMat2) ask

    case choice of
        DifferentLeft  ->
            return
                $ (rows m1 /. (rows m1 + rows m2)) * (1 /. Set.size neighbors)
        DifferentRight ->
            return
                $ (rows m2 /. (rows m1 + rows m2)) * (1 /. Set.size neighbors)

-- | Update the transition probability in the matrix.
updateTransProbMat :: [((Int, Int), Double)] -> WalkerCSRW ()
updateTransProbMat updates = do
    let allUpdates = updates ++ fmap (over _1 swap) updates
    oldState <- get
    put
        . WalkerStateCSRW
        . over _3 ( TransProbMatrix
                  . (\x -> accum x const allUpdates)
                  . unTransProbMatrix
                  )
        . unWalkerStateCSRW
        $ oldState

-- | Update the current position of the walker.
updateCs :: Int -> Int -> WalkerCSRW ()
updateCs c1 c2 = do
    oldState <- get
    put
        . WalkerStateCSRW
        . set _2 c2
        . set _1 c1
        . unWalkerStateCSRW
        $ oldState

-- | Choose which network to traverse and update the matrix.
diffUpdateTransProbMat :: Set.Set Int
                       -> Set.Set Int
                       -> WalkerCSRW ()
diffUpdateTransProbMat n1 n2 = do
    c1 <- fmap (view _1 . unWalkerStateCSRW) get
    c2 <- fmap (view _2 . unWalkerStateCSRW) get

    m1 <- fmap edgeMat1 ask
    m2 <- fmap edgeMat2 ask

    whichNetwork <- liftIO $ R.sample (R.Uniform 0 1 :: R.Uniform Int)

    if whichNetwork == (0 :: Int)
        then do
            x  <- liftIO . R.sample . R.randomElement . Set.toList $ n1
            ps <- sequence
                . fmap (getTransProb DifferentLeft . flip getNeighbors m1)
                . Set.toList
                $ n1
            updateTransProbMat . zip (repeat (x, c2)) $ ps
            updateCs x c2
        else do
            y  <- liftIO . R.sample . R.randomElement . Set.toList $ n2
            ps <- sequence
                . fmap (getTransProb DifferentRight . flip getNeighbors m2)
                . Set.toList
                $ n2
            updateTransProbMat . zip (repeat (c1, y)) $ ps
            updateCs c1 y

-- | Reset the walker state by choosing a random ID.
randomJump :: WalkerCSRW ()
randomJump = do
    size <- fmap (rows . unEdgeSimMatrix . edgeMat1) ask

    newC <- liftIO $ R.sample (R.Uniform 0 (size - 1) :: R.Uniform Int)
    oldState <- get
    put
        . WalkerStateCSRW
        . set _2 newC
        . set _1 newC
        . unWalkerStateCSRW
        $ oldState

-- | Get the transition probability matrix from a pair of similarity
-- matrices.
getTransProbMat :: Counter -> Counter -> WalkerCSRW ()
getTransProbMat counterStop counter = do
    m1 <- fmap edgeMat1 ask
    m2 <- fmap edgeMat2 ask

    c1 <- fmap (view _1 . unWalkerStateCSRW) get
    c2 <- fmap (view _2 . unWalkerStateCSRW) get

    let n1      = getNeighbors c1 m1
        n2      = getNeighbors c2 m2
        overlap = Set.union n1 n2

    if c1 == c2
        then do
            newC <- liftIO . R.sample . R.randomElement . Set.toList $ overlap
            ps   <- sequence
                  . fmap (const (getTransProb Same overlap))
                  . Set.toList
                  $ overlap
            updateTransProbMat
                . zip (zip (Set.toList overlap) (Set.toList overlap))
                $ ps
            updateCs newC newC
        else do
            diffUpdateTransProbMat n1 n2

    jump    <- liftIO $ R.sample (R.Uniform 0 1 :: R.Uniform Double)
    restart <-
        fmap
            (unWalkerRestart . (restart :: EnvironmentCSRW -> WalkerRestart))
            ask

    if jump <= restart
        then randomJump
        else return ()

    unless (counter > counterStop)
        . getTransProbMat counterStop
        $ (counter + Counter 1)

-- | Get the node correspondence scores from a transition probability
-- matrix. simVec is normalized to sum to 1.
getWalkerNodeCorrespondenceScores :: WalkerRestart
                                  -> SimVector
                                  -> TransProbMatrix
                                  -> NodeCorrScores
getWalkerNodeCorrespondenceScores (WalkerRestart restart) (SimVector simVec) =
    NodeCorrScores
        . VS.convert
        . (+ (scalar restart * (VS.map (/ VS.sum simVec) simVec)))
        . (* scalar (1 - restart))
        . largestLeftEig
        . eig
        . unTransProbMatrix

-- | Run the random walker.
evalWalker :: VertexSimMap
           -> EdgeSimMap
           -> WalkerRestart
           -> Counter
           -> LevelName
           -> LevelName
           -> IO NodeCorrScores
evalWalker vMap (EdgeSimMap eMap) restart counterStop l1 l2 = do
    let randomC = R.Uniform 0 (size - 1) :: R.Uniform Int

    c1 <- R.sample randomC
    c2 <- R.sample randomC

    transProbMat <- fmap (view _3 . unWalkerStateCSRW)
                  . flip execStateT (st c1 c2)
                  . runReaderT
                        (unWalkerCSRW $ getTransProbMat counterStop (Counter 0))
                  $ env

    let scores = getWalkerNodeCorrespondenceScores
                    restart
                    (SimVector . takeDiag . unVertexSimMatrix $ vertexSim)
                    transProbMat

    return scores
  where
    st !c1 !c2 =
        WalkerStateCSRW (c1, c2, TransProbMatrix . (size><size) . repeat $ 0)
    size       = rows . unEdgeSimMatrix . edgeMat1 $ env
    env = EnvironmentCSRW { edgeMat1 = lookupLevel l1 eMap
                          , edgeMat2 = lookupLevel l2 eMap
                          , vertexMat = vertexSim
                          , restart = restart
                          }
    vertexSim     = getVertexSim l1 l2 vMap
    lookupLevel l = lookupWithError
                        (levelErr l)
                        l
    levelErr l    = ("Level: " ++ (T.unpack $ unLevelName l) ++ " not found.")

