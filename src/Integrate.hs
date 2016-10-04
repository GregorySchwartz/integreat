{- Integrate
Gregory W. Schwartz

Collections the functions pertaining to the integration of levels of data.
-}

{-# LANGUAGE BangPatterns #-}

module Integrate
    ( integrateCosineSim
    , integrateWalker
    , integrateCSRW
    ) where

-- Standard
import qualified Data.Set as Set
import qualified Data.Map.Strict as Map

-- Cabal
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Data.Text as T
import Data.Graph.Inductive

-- Local
import Types
import Utility
import Alignment.Cosine
import Alignment.RandomWalk
import Alignment.CSRW

-- | Do integration using cosine similarity.
integrateCosineSim :: VertexSimMap
                   -> EdgeSimMap
                   -> IO NodeCorrScores
integrateCosineSim vertexSimMap (EdgeSimMap edgeSimMap) =
    return
        . NodeCorrScores
        . V.convert
        . applyRows avgVec
        . fmap unNodeCorrScores
        $ compare <$> levels <*> levels
  where
    compare x y   =
        cosineIntegrate vertexSimMap x y (lookupLevel x) (lookupLevel y)
    levels        = Map.keys edgeSimMap
    lookupLevel l = lookupWithError
                        ( "Level: "
                       ++ (T.unpack $ unLevelName l)
                       ++ " notFound in integration step."
                        )
                        l
                        edgeSimMap

-- | Do integration using a random walker.
integrateWalker :: WalkerRestart -> Counter -> GrMap -> IO NodeCorrScores
integrateWalker walkerRestart counterStop (GrMap grMap) = do
    fmap ( NodeCorrScores
         . V.convert
         . applyRows avgVec
         . fmap VS.fromList
         )
        . sequence
        $ eval <$> Map.elems grMap <*> Map.elems grMap
  where
    eval gr1 gr2 = sequence
                 . fmap (randomWalkerScore walkerRestart counterStop gr1 gr2)
                 . allNodes (unLevelGr gr1)
                 . unLevelGr
                 $ gr2
    allNodes gr1 gr2 = Set.toList
                     . Set.union (Set.fromList . nodes $ gr2)
                     . Set.fromList
                     . nodes
                     $ gr1

-- | Use integration using a random walker "CSRW".
integrateCSRW
    :: VertexSimMap
    -> EdgeSimMap
    -> WalkerRestart
    -> Counter
    -> IO NodeCorrScores
integrateCSRW vertexSimMap edgeSimMap restart counterStop = do
      fmap ( NodeCorrScores
           . V.convert
           . applyRows avgVec
           . fmap unNodeCorrScores
           )
        . sequence
        $ eval <$> levels <*> levels
  where
    eval           = evalWalker vertexSimMap edgeSimMap restart counterStop
    levels         = Map.keys . unEdgeSimMap $ edgeSimMap
    lookupLevel l  = lookupWithError
                        ("Level: " ++ (T.unpack $ unLevelName l) ++ " notFound")
                        l
                   . unEdgeSimMap
                   $ edgeSimMap
