
class BuildingBlock:

    ### DUNDERS

    def __init__(self, smiles, molecular_weight=None, catalog_metadata=None):
        
        self._smiles = smiles
        self._molecular_weight = molecular_weight
        self._catalog_metadata = catalog_metadata

    def __repr__(self):
        return f'BuildingBlock({self.smiles})'

class BuildingBlockSet(set):

    pass

