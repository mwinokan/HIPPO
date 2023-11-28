
import mout
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
        self._chain = None 
        self._site_index = None 
        self._pose_name = None
        self._crystal_name = None
        self._alternate_name = None
        self._site_name = None
        self._pdb_entry = None

        # arguments
        self.name = name
        self.smiles = smiles
        self._tags = tags or []

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

        self.__init__(name, metadata['new_smiles'] or metadata['smiles'])

        self.chain = chain
        self.site_index = site_index

        self.pose_name = metadata['crystal_name']
        self.crystal_name = metadata['RealCrystalName']
        self.alternate_name = metadata['alternate_name']
        self.site_name = metadata['site_name']
        self.pdb_entry = metadata['pdb_entry']

        ### CREATE THE POSE OBJECT

        pose = Pose.from_bound_pdb(self, path, metadata)
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
            mout.error(f'There is a dot in the SMILES [{self.name}]: {s}')
            s = sorted(s.split('.'), key=lambda x: len(x))[-1]
            mout.warning(f'Using: {s}')
        
        # flatten the smiles
        self._stereo_smiles = s
        self._smiles = s.replace('@','')

    @property
    def chain(self):
        return self._chain

    @chain.setter
    def chain(self,a):
        assert isinstance(a,str)
        assert len(a) == 1
        self._chain = a

    @property
    def site_index(self):
        return self._site_index

    @site_index.setter
    def site_index(self,a):
        a = int(a)
        self._site_index = a
    
    @property
    def pose_name(self):
        return self._pose_name

    @pose_name.setter
    def pose_name(self,a):
        assert isinstance(a,str)
        self._pose_name = a

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
    def site_name(self):
        return self._site_name

    @site_name.setter
    def site_name(self,a):
        assert isinstance(a,str)
        self._site_name = a

    @property
    def pdb_entry(self):
        return self._pdb_entry

    @pdb_entry.setter
    def pdb_entry(self,a):
        # assert isinstance(a,str)
        self._pdb_entry = a

    @property
    def poses(self):
        return self._poses

    @property
    def num_poses(self):
        return len(self.poses)
    
### METHODS

    def add_pose(self, pose):
        self.poses.append(pose)

    def summary(self):

        mout.header(str(self))

        mout.var('chain',self.chain)
        mout.var('site_index',self.site_index)
        mout.var('pose_name',self.pose_name)
        mout.var('crystal_name',self.crystal_name)
        mout.var('pdb_entry',self.pdb_entry)
        mout.var('poses',[p.name for p in self.poses])

### DUNDERS

    def __repr__(self):
        return f'Compound({self.name}, {self.smiles}, #poses={self.num_poses})'

    def __eq__(self, other):
        return self.name == other.name
