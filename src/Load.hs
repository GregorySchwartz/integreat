{- Load
Gregory W. Schwartz

Collections the functions pertaining to the loading and unification of data
sets.
-}

{-# LANGUAGE BangPatterns #-}

module Load
    ( datasToEntities
    , entitiesToLevels
    , unifyAllLevels
    , filterEntities
    , getIDVec
    , getIDMap
    , defVertexSimMap
    , vertexCsvToLevels
    , standardizeLevel
    , getPremadeNetworks
    ) where

-- Standard
import Data.List
import qualified Data.Sequence as Seq
import qualified Data.Set as Set
import qualified Data.Map.Strict as Map
import qualified Data.Foldable as F
import Data.Function (on)

-- Cabal
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Data.Text as T
import Data.Graph.Inductive
import qualified Control.Lens as L
import Numeric.LinearAlgebra

-- Local
import Types
import Utility

-- | Convert data entries to entities.
datasToEntities :: [DataEntry] -> [Entity]
datasToEntities = fmap dataEntryToEntity

-- | Convert a vertex entry to an entity.
dataEntryToEntity :: DataEntry -> Entity
dataEntryToEntity x = Entity { _entityID    = ID $ vertex x
                             , _dataSetName = DataSetName $ dataReplicate x
                             , _levelName   = LevelName $ dataLevel x
                             , _entityValue = intensity x
                             }

-- | Convert a collection of entities to levels.
entitiesToLevels :: [Entity] -> [(LevelName, Level)]
entitiesToLevels = fmap (L.over L._2 Level)
                 . Map.toAscList
                 . Map.unionsWith (Map.unionWith Map.union)
                 . fmap (\ !x -> Map.singleton
                                    (_levelName x)
                                    ( Map.singleton
                                        (_entityID x)
                                        (Map.singleton (_dataSetName x) x)
                                    )
                        )

-- | Combine all levels together to get a unique list of IDs.
unifyAllLevels :: [Level] -> UnifiedData
unifyAllLevels = UnifiedData . Map.unionsWith Map.union . fmap unLevel

-- | Filter entities that appear in less than the specified amount.
filterEntities :: NumSamples -> [Entity] -> [Entity]
filterEntities (NumSamples numSamples) = mconcat
                                       . filter ((>= numSamples) . length)
                                       . groupBy ((==) `on` _entityID)
                                       . sortBy (compare `on` _entityID)
    
-- | Get a vector of all IDs for easy indexing.
getIDVec :: UnifiedData -> IDVec
getIDVec =
    IDVec . (\x -> V.fromList . Map.keys $ x) . unUnifiedData

-- | Get a map of all IDs for conversion.
getIDMap :: UnifiedData -> IDMap
getIDMap = IDMap
         . Map.fromList
         . V.toList
         . V.imap (\i v -> (v, i))
         . unIDVec
         . getIDVec

-- | The default vertex similarity map using just a diagonal of 1 matrix
-- (same entities are 1, different are 0).
defVertexSimMap :: IDMap -> [LevelName] -> VertexSimMap
defVertexSimMap (IDMap idMap) xs = VertexSimMap
                                 . Map.fromListWith Map.union
                                 . fmap ( \((!x, !y), !z)
                                       -> (x, Map.singleton y z)
                                        )
                                 . flip zip (repeat defSimMat)
                                 $ (,) <$> xs <*> xs
  where
    defSimMat = VertexSimMatrix
              . diag
              . flip VS.replicate (1 :: Double)
              . Map.size
              $ idMap

-- | Convert a list of vertex entries to a vertex similarity map by joining
-- together all like comparisons and converting to those values into
-- matrices.
vertexCsvToLevels :: IDMap -> [VertexEntry] -> VertexSimMap
vertexCsvToLevels (IDMap idMap) =
    VertexSimMap
        . (Map.map . Map.map) ( VertexSimMatrix
                              . assoc (Map.size idMap, Map.size idMap) 0
                              . F.toList
                              )
        . Map.unionsWith (Map.unionWith (Seq.><))
        . fmap toLevelSimMap
  where
    toLevelSimMap x = Map.singleton
                        (LevelName $ vertexLevel1 x)
                        ( Map.singleton (LevelName $ vertexLevel2 x)
                                        ( Seq.singleton ( ( idLookup $ vertex1 x
                                                          , idLookup $ vertex2 x
                                                          )
                                                        , similarity x
                                                        )
                                        )
                        )
    idLookup x =
        lookupWithError ("ID: " ++ T.unpack x ++ " not found.") (ID x) idMap

-- | Convert a level's IDs to be integers based on the universal labeling.
standardizeLevel :: IDMap -> Level -> StandardLevel
standardizeLevel (IDMap idMap) (Level level) =
    StandardLevel
        . Map.map (flip standardizeDataSets standardDataSets)
        . Map.mapKeys (\ !x -> (x, lookupWithError (keyNotFound x) x idMap))
        $ level
  where
    standardDataSets = getStandardDataSets . Level $ level
    keyNotFound k    = "ID: " ++ show k ++ " not found."

-- | Unify the ordering of the data for data set entries.
getStandardDataSets :: Level -> StandardDataSets
getStandardDataSets = StandardDataSets
                    . Seq.fromList
                    . Set.toList
                    . Set.fromList
                    . concatMap Map.keys
                    . Map.elems
                    . unLevel

-- | Unify the ordering of the data for data set entries.
standardizeDataSets :: Map.Map DataSetName Entity
                    -> StandardDataSets
                    -> (Seq.Seq (Maybe Entity))
standardizeDataSets m = fmap (flip Map.lookup m) . unStandardDataSets

-- | Get a vector of all IDs from pre-made networks.
getPremadeIDVec :: [String] -> IDVec
getPremadeIDVec = IDVec . V.fromList . fmap (ID . T.pack)

-- | Get a map of all IDs for conversion. Assumes the IDs are already ordered
-- correctly.
getPremadeIDMap :: [String] -> IDMap
getPremadeIDMap xs =
    IDMap . Map.fromList . flip zip (fmap read xs) . fmap (ID . T.pack) $ xs

-- | Get the graph from a list of edges. Assumes the IDs integers,
-- starting at 0.
getPremadeGr :: IDMap -> [(String, String, Double)] -> LevelGr
getPremadeGr (IDMap idMap) =
    LevelGr
        . undir
        . mkGraph (zip [0..Map.size idMap] [0..Map.size idMap])
        . fmap (\(!x, !y, !z) -> (read x, read y, z))

-- | Get the edge matrix from a list of edges. Assumes the IDs integers,
-- starting at 0.
getPremadeEdgeSimMatrix :: IDMap -> [(String, String, Double)] -> EdgeSimMatrix
getPremadeEdgeSimMatrix (IDMap idMap) =
    EdgeSimMatrix
        . assoc (Map.size idMap, Map.size idMap) 0
        . fmap (\(!x, !y, !z) -> ((read x, read y), z))

-- | Get the pre-made level networks for use with integration. We assume that
-- the IDs are integers, starting at 0, with everything in order.
getPremadeNetworks
    :: (Maybe [String], [([String], [(String, String, Double)])])
    -> (Maybe (Set.Set ID), Maybe UnifiedData, IDMap, IDVec, VertexSimMap, EdgeSimMap, GrMap)
getPremadeNetworks (!changedVerticesString, !allNetworks) =
    (changedVertices, Nothing, idMap, idVec, vertexSimMap, edgeSimMap, grMap)
  where
    changedVertices =
        fmap (Set.fromList . fmap (ID . T.pack)) $ changedVerticesString
    grMap        = GrMap
                 . Map.fromList
                 . zip levelNames
                 . fmap (getPremadeGr idMap)
                 . fmap snd
                 $ allNetworks
    edgeSimMap   = EdgeSimMap
                 . Map.fromList
                 . zip levelNames
                 . fmap (getPremadeEdgeSimMatrix idMap)
                 . fmap snd
                 $ allNetworks
    vertexSimMap = defVertexSimMap idMap levelNames
    levelNames   = fmap (LevelName . T.pack . show) [0..length allNetworks - 1]
    idMap        = getPremadeIDMap allVertices
    idVec        = getPremadeIDVec allVertices
    allVertices  = concatMap fst allNetworks
