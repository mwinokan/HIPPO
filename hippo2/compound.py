
import mout
import mcol
import molparse as mp
from rdkit import Chem
from pathlib import Path

from .pose import Pose

class Compound:

    def __init__(self, name, smiles, tags=None):

        # blanks
        self._name = None
        self._smiles = None
        self._stereo_smiles = None
        self._orig_smiles = None
        self._mol = None
        self._poses = []

        # from XChem crystal structure metadata
        self._crystal_name = None
        self._alternate_name = None

        # arguments
        self.name = name
        self.smiles = smiles
        self.tags = list(tags or [])

### FACTORIES

    @classmethod
    def from_mol(cls, name, path, tags=None):

        assert isinstance(path, Path)
    
        tags = tags or []
        if isinstance(tags,str):
            tags = [tags]
        assert isinstance(tags,list)

        mol = Chem.MolFromMolFile(path)
        smiles = mp.mol_to_smiles(mol)
        
        self = cls.__init__(name, smiles, tags)

        return self

    @classmethod
    def from_bound_pdb(cls, name, path, metadata, site_index=None, chain=None, tags=None):

        assert isinstance(path, Path)

        pose_name = metadata['crystal_name']

        if site_index is None:
            site_index = pose_name[-2]

        if chain is None:
            chain = pose_name[-1]

        tags = tags or []   

        self = cls.__new__(cls)

        self.__init__(name, metadata['new_smiles'] or metadata['smiles'], tags)
        
        self.crystal_name = metadata['RealCrystalName']
        self.alternate_name = metadata['alternate_name']
        
        ### CREATE THE POSE OBJECT

        pose = Pose.from_bound_pdb(self, path, metadata, tags=tags)
        self.add_pose(pose)

        return self

### PROPERTIES

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self,a):
        assert self._name is None
        self._name = a

    @property
    def smiles(self):
        return self._smiles

    @smiles.setter
    def smiles(self,s):

        if self._smiles is not None:
            mout.warning(f'Overwriting SMILES: {self.name}')
        assert isinstance(s,str)

        self._orig_smiles = s

        # if multiple molecules take the largest
        if '.' in s:
            s = sorted(s.split('.'), key=lambda x: len(x))[-1]
        
        # flatten the smiles
        self._stereo_smiles = s
        self._smiles = s.replace('@','')

        if self._smiles != self._orig_smiles:

            annotated_smiles_str = self.orig_smiles.replace('.',f'{mcol.error}{mcol.underline}.{mcol.clear}{mcol.warning}')
            annotated_smiles_str = annotated_smiles_str.replace('@',f'{mcol.error}{mcol.underline}@{mcol.clear}{mcol.warning}')

            mout.warning(f'SMILES was changed: {annotated_smiles_str} --> {self.smiles}')

    @property
    def orig_smiles(self):
        return self._orig_smiles

    @property
    def orig_smiles(self):
        return self._orig_smiles

    @property
    def crystal_name(self):
        return self._crystal_name

    @crystal_name.setter
    def crystal_name(self,a):
        assert isinstance(a,str)
        self._crystal_name = a

    @property
    def alternate_name(self):
        return self._alternate_name

    @alternate_name.setter
    def alternate_name(self,a):
        assert isinstance(a,str)
        self._alternate_name = a

    @property
    def poses(self):
        return self._poses

    @property
    def num_poses(self):
        return len(self.poses)

    @property
    def mol(self):
        if self._mol is None:
            self._mol = Chem.MolFromSmiles(self.smiles)
        return self._mol
    
### METHODS

    def add_pose(self, pose):
        self.poses.append(pose)

    def summary(self):

        mout.header(str(self))

        mout.var('crystal_name',self.crystal_name)
        mout.var('poses',[p.name for p in self.poses])

### DUNDERS

    def __repr__(self):
        return f'Compound({self.name}, {self.smiles}, #poses={self.num_poses})'

    def __eq__(self, other):
        return self.name == other.name
