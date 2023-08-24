
import os
import mout
import molparse as mp
from .compound import Compound
from .tools import df_row_to_dict
import pandas as pd

class CompoundSet:

    ### DUNDERS

    def __init__(self, name, compounds = None):
        
        self._name = name
        if compounds is None:
            compounds = []
        self._compounds = compounds
        self._fingerprint_df = None
        self._metadata_cols = []

    def __repr__(self):
        return f'CompoundSet("{self.name}", #compounds={self.num_compounds})'

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
    
    @property
    def compounds(self):
        return self._compounds
   
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
            self._fingerprint_df = pd.DataFrame([f for f in self.fingerprints])
        return self._fingerprint_df
    
    ### METHODS

    def add_compound(self, compound):
        self.compounds.append(compound)
        compound._set_name = self.name

    def get_present_features(self):

        features = set()
        for col in [c for c in self.fingerprint_df.columns if c not in self._metadata_cols]:
            covered = any(self.fingerprint_df[col].values)
            if covered:
                features.add(col)

        return features
