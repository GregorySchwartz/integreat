{- Integrate
Gregory W. Schwartz

Collections the functions pertaining to the integration of levels of data.
-}

{-# LANGUAGE BangPatterns #-}

module Integrate
    ( integrate
    ) where

-- Standard
import qualified Data.Map.Strict as Map

-- Cabal
import qualified Data.Vector as V
import qualified Data.Text as T

-- Local
import Types
import Utility
import Alignment.Cosine
import Alignment.RandomWalk

-- | Get the scores of each vertex for the best conserved vertices (based on the
-- vertex and edge similarity in each level network).
integrate :: AlignmentMethod
          -> VertexSimMap
          -> EdgeSimMap
          -> WalkerRestart
          -> Counter
          -> IO NodeCorrScores
integrate CosineSimilarity vertexSimMap (EdgeSimMap edgeSimMap) _ _ =
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
integrate RandomWalker vertexSimMap edgeSimMap restart counterStop = do
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
