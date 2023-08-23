
import os
import molparse as mp
from .compound import Compound
from .tools import df_row_to_dict

class CompoundSet:

    ### DUNDERS

    def __init__(self, name, compounds):
        
        self._name = name
        self._compounds = compounds

    def __repr__(self):
        return f'CompoundSet("{self.name}", #compounds={self.num_compounds})'

    ### FACTORIES

    @classmethod
    def from_df(cls, name, df):

        self = cls.__new__(cls)

        id_col = df.columns[0]
        mol_col = df.columns[1]

        compounds = []
        for index, row in df.iterrows():
            mol_name = row[id_col]
            mol = row[mol_col]
            compound = Compound.from_rdkit_mol(mol_name,mol)
            compound._set_name = name
            compounds.append(compound)
        
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
    
