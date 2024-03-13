
import os
import mout
import molparse as mp
from .compound import Compound
from .tools import df_row_to_dict
import pandas as pd
import random

from collections.abc import MutableSet

class CompoundSet(MutableSet):

    ### DUNDERS

    def __init__(self, name, compounds=(), immutable=False):
        # super(CompoundSet,self).__init__(compounds)

        self._name = name
        self._elements = []
        self._immutable = immutable

        for comp in compounds:
            self.add(comp)

        self._fingerprint_df = None
        self._metadata_cols = []

    def __repr__(self):
        bbs = self.get_building_blocks()
        # purchaseable_bbs = self.get_building_blocks(purchaseable_only=True)
        # purchaseable_bbs.get_products(self.compounds)
        # return f'CompoundSet("{self.name}", #compounds={self.num_compounds}, #bbs={len(bbs)}, #bbs (purchaseable)={len(purchaseable_bbs)}: ${purchaseable_bbs.get_price():.2f})'
        return f'CompoundSet("{self.name}", #compounds={self.num_compounds}, #bbs={len(bbs)})'

    def __iter__(self):
        return iter(self._elements)

    def __contains__(self, compound):
        return compound in self._elements

    def __len__(self):
        return len(self._elements)

    def __getitem__(self, key):

        if isinstance(key,list):
            return CompoundSet('queried',[c for c in [self[k] for k in key] if c is not None])

        if isinstance(key,Compound):
            key = str(key)

        elif isinstance(key,int):
            return [c for c in self][key]

        matches = [comp for comp in self if str(comp) == key]

        if len(matches) < 1:
            mout.error(f'{key} not in {self}')
            return None
        elif len(matches) > 1:
            mout.error(f'Multiple {key} in {self}')
            exit()
            return None

        return matches[0]

    def __setitem__(self, key, value):
        assert isinstance(key,int)
        [c for c in self][key] = value

    ### FACTORIES

    @classmethod
    def from_df(cls, name, df, protein, verbosity=1):

        self = cls.__new__(cls)

        id_col = df.columns[0]
        mol_col = df.columns[1]

        compounds = []
        for index, row in df.iterrows():
            mol_name = row[id_col]
            if verbosity:
                mout.progress(index,len(df),prepend='mol --> Compound',append=mol_name, fill='🦛')
            mol = row[mol_col]
            compound = Compound.from_rdkit_mol(mol_name,mol)
            compound._set_name = name
            compound._protein_system = protein
            compounds.append(compound)

        if verbosity:
            mout.finish()
        
        self.__init__(name,compounds)

        return self

    @classmethod
    def from_bound_pdbs(cls, name, pdbs, metadata_df, prefix):

        self = cls.__new__(cls)

        compounds = []
        for pdb in pdbs:
            pose_name = os.path.basename(pdb).removesuffix('_bound.pdb')
            comp_name = pose_name.removeprefix(prefix)
            metadata = df_row_to_dict(metadata_df[ metadata_df['crystal_name'] == pose_name ])
            compound = Compound.from_bound_pdb(comp_name, pdb, metadata)
            compound._set_name = name
            compounds.append(compound)

        self.__init__(name,compounds)

        return self

    ### PROPERTIES

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self,a):
        self._name = a
    
    @property
    def compounds(self):
        return self._elements
   
    @property
    def num_compounds(self):
        return len(self.compounds)

    @property
    def fingerprinted_compounds(self):
        return [c for c in self.compounds if c.fingerprint is not None]

    @property
    def fingerprints(self):
        return [c.fingerprint for c in self.compounds if c.fingerprint is not None]

    @property
    def fingerprint_df(self):
        if self._fingerprint_df is None:
            fingerprints = self.fingerprints
            if len(fingerprints) < 1:
                mout.error(f'no fingerprints for {self}')
            self._fingerprint_df = pd.DataFrame(fingerprints)
        return self._fingerprint_df
    
    @property
    def immutable(self):
        return self._immutable
    
    @immutable.setter
    def immutable(self,b):
        self._immutable = b

    ### METHODS

    def pop(self):
        return self._elements.pop()
        # last = self._elements[-1]
        # self._elements = self._elements[:-1]
        # return last

    def discard(self, key):
        assert not self.immutable
        if key in self:
            i = self._elements.index(key)
            del self._elements[i]
        else:
            raise ValueError(f'{key} not in {self}')

    def remove(self, key):
        assert not self.immutable
        if key in self:
            i = self._elements.index(key)
            del self._elements[i]
        else:
            raise ValueError(f'{key} not in {self}')

    def add(self, compound):
        assert not self.immutable
        if compound not in self._elements:
            self._elements.append(compound)
            compound._set_name = self.name
        else:
            raise ValueError(f'{compound} already in {self}')

    def shuffle(self):
        random.shuffle(self._elements)

    # def get_random(size=1):
    #     return random.sample(self._elements,size)

    def get_present_features(self):

        features = set()
        for col in [c for c in self.fingerprint_df.columns if c not in self._metadata_cols]:
            covered = any(self.fingerprint_df[col].values)
            if covered:
                features.add(col)

        return features

    def summary(self):

        mout.header(self)

        print_data = []
        for comp in self:
            print_data.append(dict(
                smiles=comp.smiles,
                set_name=comp.set_name,
                is_pains=comp.is_pains,
                building_blocks=comp.building_blocks,
                cost_range_str=comp.cost_range_str,
                lead_time=comp.lead_time,
            ))

        print(pd.DataFrame(print_data))

    def get_building_blocks(self,purchaseable_only=False):

        from .block import BuildingBlockSet

        bb_set = BuildingBlockSet()

        for comp in self:
            if comp.building_blocks is not None:
                for bb in comp.building_blocks:
                    if not purchaseable_only or bb.purchaseable:
                        bb_set.add(bb)

        return bb_set

    ### PROTECTION

    def __le__(self):
        raise NotImplementedError('CompoundSet.__le__')

    def __lt__(self):
        raise NotImplementedError('CompoundSet.__lt__')

    def __eq__(self):
        raise NotImplementedError('CompoundSet.__eq__')

    def __ne__(self):
        raise NotImplementedError('CompoundSet.__ne__')

    def __gt__(self):
        raise NotImplementedError('CompoundSet.__gt__')

    def __ge__(self):
        raise NotImplementedError('CompoundSet.__ge__')

    def __and__(self):
        raise NotImplementedError('CompoundSet.__and__')

    def __or__(self):
        raise NotImplementedError('CompoundSet.__or__')

    # def __sub__(self):
    #     raise NotImplementedError('CompoundSet.__sub__')

    def __xor__(self):
        raise NotImplementedError('CompoundSet.__xor__')

    def isdisjoint(self):
        raise NotImplementedError('CompoundSet.isdisjoint')

    def __ior__(self):
        raise NotImplementedError('CompoundSet.__ior__')

    def __ixor__(self):
        raise NotImplementedError('CompoundSet.__ixor__')

    def __iand__(self):
        raise NotImplementedError('CompoundSet.__iand__')

    def __isub__(self):
        raise NotImplementedError('CompoundSet.__isub__')
