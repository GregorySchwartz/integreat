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
    , getIDVec
    , getIDMap
    , defVertexSimMap
    , vertexCsvToLevels
    , standardizeLevel
    ) where

-- Standard
import qualified Data.Sequence as Seq
import qualified Data.Map.Strict as Map
import qualified Data.Foldable as F

-- Cabal
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Data.Text as T
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

-- | Convert an collection of entities to levels.
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
standardizeLevel (IDMap idMap) =
    StandardLevel
        . Map.mapKeys (\ !x -> (x, lookupWithError (keyNotFound x) x idMap))
        . unLevel
  where
    keyNotFound k = "ID: " ++ show k ++ " not found."
