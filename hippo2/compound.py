
import mout
import mcol
import molparse as mp
from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem.rdMolHash import MolHash, HashFunction
from pathlib import Path

from .pose import Pose
from .tset import TagSet
from .pset import PoseSet
from .rset import ReactionSet
from .tools import remove_isotopes_from_smiles

import re

class Compound:

    def __init__(self, name, smiles, tags=None):

        # mout.debug(f'Compound.__init__({name})')

        # blanks
        self._name = None
        self._smiles = None
        self._stereo_smiles = None
        self._orig_smiles = None
        self._mol = None
        self._poses = PoseSet()
        self._inspirations = []
        self._base = None
        self._fp_1024 = None
        self._canonical_smiles_hash = None
        # self._amount = None # mg
        self._reactions = ReactionSet()

        # from XChem crystal structure metadata
        self._crystal_name = None
        self._alternate_name = None

        # arguments
        self.name = name
        self.smiles = smiles
        self.tags = TagSet(tags or [])

### FACTORIES

    @classmethod
    def from_mol(cls, name, path, tags=None):

        if isinstance(path, Path):
            mol = Chem.MolFromMolFile(str(path))
        else:
            mol = path
    
        tags = tags or []
        if isinstance(tags,str):
            tags = [tags]
        tags = TagSet(tags)

        smiles = mp.rdkit.mol_to_smiles(mol)
        
        self = cls.__new__(cls)

        self.__init__(name, smiles, tags)

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
        tags = TagSet(tags)

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

        # remove isotopic stuff
        if re.search(r'([\[][0-9]+[A-Z]+\])',self._smiles):
            mout.warning(f'Isotope(s) in SMILES: {self._smiles}')
            self._smiles = remove_isotopes_from_smiles(self._smiles)

        return

        if self._smiles != self._orig_smiles:

            annotated_smiles_str = self.orig_smiles.replace('.',f'{mcol.error}{mcol.underline}.{mcol.clear}{mcol.warning}')
            annotated_smiles_str = annotated_smiles_str.replace('@',f'{mcol.error}{mcol.underline}@{mcol.clear}{mcol.warning}')
            
            # if re.search(r'([\[][0-9]+[A-Z]+\])',self._smiles):
            #     print(self._smiles)

            #     for match in set(re.findall(r'([\[][0-9]+[A-Z]+\])',self._smiles)):
            #         print(match)
            #         annotated_smiles_str.replace(match, f'{mcol.error}{mcol.underline}{match}{mcol.clear}{mcol.warning}')

            mout.warning(f'SMILES was changed: {annotated_smiles_str} --> {self.smiles}')

    @property
    def orig_smiles(self):
        return self._orig_smiles

    @property
    def stereo_smiles(self):
        return self._stereo_smiles

    @property
    def base(self):
        return self._base

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
    def pose(self):
        if self.num_poses > 1:
            mout.warning('Returning first of multiple poses')
        return self.poses[0]

    @property
    def reaction(self):
        assert self.num_reactions == 1
        return self.reactions[0]

    @property
    def reactions(self):
        return self._reactions

    @property
    def num_reactions(self):
        return len(self.reactions)

    @property
    def building_blocks(self):
        return self.reaction.reactants

    @property
    def mol(self):
        if self._mol is None:
            self._mol = Chem.MolFromSmiles(self.smiles)
        return self._mol

    @property
    def inspirations(self):
        return self._inspirations
    
    @inspirations.setter
    def inspirations(self, s):
        self._inspirations = s

    @property
    def fp_1024(self):
        if self._fp_1024 is None:
            self._fp_1024 = FingerprintMols.FingerprintMol(self.mol, minPath=1, maxPath=7, fpSize=2048, bitsPerHash=2, useHs=True, tgtDensity=0.0, minSize=128)
        return self._fp_1024

    @property
    def canonical_smiles_hash(self):
        if self._canonical_smiles_hash is None:
            self._canonical_smiles_hash = hash(MolHash(self.mol, HashFunction.CanonicalSmiles))
        return self._canonical_smiles_hash
    
    # @property
    # def dict(self):
    #     return dict(
    #         name = self.name,
    #         smiles = self.smiles,
    #         orig_smiles = self.orig_smiles,
    #         crystal_name = self.crystal_name,
    #         crystal_name = self.crystal_name,

    #     )

### METHODS

    def add_pose(self, pose):
        self.poses.add(pose)

    def add_reaction(self, reaction):
        self.reactions.add(reaction)

    def summary(self):

        mout.header(str(self))

        mout.var('smiles',self.smiles)
        
        if self.orig_smiles: mout.var('orig_smiles',self.orig_smiles)
        if self.stereo_smiles: mout.var('stereo_smiles',self.stereo_smiles)
        if self.crystal_name: mout.var('crystal_name',self.crystal_name)
        if self.inspirations: mout.var('inspirations',self.inspirations)
        # if self.amount is not None: mout.var('amount',self.amount)
        if self.base: mout.var('base',self.base)
        if self.alternate_name: mout.var('alternate_name',self.alternate_name)
        
        mout.var('tags',self.tags)
        mout.var('poses',self.poses)

        if self.num_reactions != 1:
            mout.var('reactions',self.reactions, symbol='')
        else:
            # mout.var('reactions:',"\n")
            mout.out("")
            self.reaction.summary()

    def get_pose(self, name):
        return [p for p in self.poses if p.name == name][0]

### DUNDERS

    def __repr__(self):
        return f'Compound({self.name}, {self.smiles}, #poses={self.num_poses})'

    def __eq__(self, other):
        return self.name == other.name

    def __hash__(self):
        return self.canonical_smiles_hash
