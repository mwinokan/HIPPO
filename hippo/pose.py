
import mcol
import molparse as mp

from rdkit import Chem

import numpy as np

from .tags import TagSet

import pickle

import logging
logger = logging.getLogger('HIPPO')

from molparse.rdkit.features import FEATURE_FAMILIES, COMPLEMENTARY_FEATURES

CUTOFF_PADDING = 1.0

FEATURE_PAIR_CUTOFFS = {
	'Donor Acceptor': 3.5 + CUTOFF_PADDING,
	'Acceptor Donor': 3.5 + CUTOFF_PADDING,
	'NegIonizable PosIonizable': 4.5 + CUTOFF_PADDING,
	'PosIonizable NegIonizable': 4.5 + CUTOFF_PADDING,
	'Aromatic PosIonizable': 4.5 + CUTOFF_PADDING,
	'PosIonizable Aromatic': 4.5 + CUTOFF_PADDING,
	'Aromatic Aromatic': 6.0 + CUTOFF_PADDING,
	'Hydrophobe Hydrophobe': 4.5 + CUTOFF_PADDING,
}

class Pose:

	"""A :class:`.Pose` is a particular conformer of a :class:`.Compound` within a protein environment. A pose will have its own (stereochemical) smiles string, and must have a path to a coordinate file. Poses can have *inspirations* that can be used to trace fragment-derived scaffolds in merges and expansions.

	:class:`.Pose` objects should not be created directly. Instead use :meth:`.HIPPO.register_pose` or :meth:`.HIPPO.poses`
	"""

	def __init__(self, 
		db,
		id: int,
		name: str,
		smiles: str,
		reference: int, # another pose
		path: str,
		compound: int,
		target: str,
		mol: Chem.Mol | bytes | None,
		fingerprint: bytes,
		metadata: dict | None = None,
	):

		self._db = db
		self._id = id
		self._name = name
		self._smiles = smiles
		self._compound_id = compound
		self._target = target
		self._path = path
		self._protein_system = None
		self._table = 'pose'

		if fingerprint:
			# logger.debug('unpickling fingerprint')
			fingerprint = pickle.loads(fingerprint)		
		self._fingerprint = fingerprint

		# print(f'{self}{metadata=}')
		self._metadata = metadata
		self._tags = None
		self._reference = reference

		if isinstance(mol, bytes):
			self._mol = Chem.Mol(mol)
		else: 
			self._mol = mol
		
	### FACTORIES

	### PROPERTIES

	@property
	def db(self):
		"""Returns a pointer to the parent database"""
		return self._db

	@property
	def id(self) -> int:
		"""Returns the pose's database ID"""
		return self._id

	@property
	def name(self) -> str:
		"""Returns the pose's name"""
		return self._name

	@property
	def smiles(self):
		"""Returns the pose's smiles"""
		if not self._smiles:
			from molparse.rdkit import mol_to_smiles
			self._smiles = mol_to_smiles(self.mol)
			self.db.update(table='pose', id=self.id, key='pose_smiles', value=self._smiles)
		return self._smiles

	@property
	def target(self):
		"""Returns the pose's associated target"""
		if isinstance(self._target, int):
			self._target = self.db.get_target(id=self._target)
		return self._target

	@property
	def compound_id(self) -> int:
		"""Returns the pose's associated compound ID"""
		return self._compound_id

	@property
	def compound(self):
		"""Returns the pose's associated compound"""
		return self.get_compound()

	@property
	def path(self):
		"""Returns the pose's path"""
		return self._path

	@property
	def reference(self):
		"""Returns the pose's protein reference (another pose)"""
		if isinstance(self._reference, int):
			self._reference = self.db.get_pose(id=self._reference)
		return self._reference

	@property
	def mol(self):
		"""Returns a pose's rdkit.Chem.Mol"""
		if not self._mol and self.path:

			if self.path.endswith('.pdb'):

				# logger.reading(self.path)
				sys = mp.parse(self.path, verbosity=False)
				
				self.protein_system = sys.protein_system
				
				lig_residues = sys['rLIG']
				if len(lig_residues) > 1:
					logger.warning(f'Multiple ligands in PDB {self}')

				lig_res = lig_residues[0]
				
				if not (mol := lig_res.rdkit_mol):
					logger.error(f'[{self}] Error computing RDKit Mol from PDB={self.path}')
					
					print(lig_res.pdb_block)

					lig_res.plot3d()

					raise Exception
					return None

				self.mol = mol

			elif self.path.endswith('.mol'):

				# logger.reading(self.path)
				mol = mp.parse(self.path, verbosity=False)
				self.mol = mol

			else:

				raise NotImplementedError

			if not mol:
				logger.error(f'Could not parse {self}.path={self.path}')

		return self._mol

	@mol.setter
	def mol(self, m):
		"""Set the pose's rdkit.Chem.Mol"""
		self._mol = m
		self.db.update(table='pose', id=self.id, key='pose_mol', value=m.ToBinary())

	@property
	def protein_system(self):
		"""Returns the pose's protein molparse.System"""
		if self._protein_system is None and self.path.endswith('.pdb'):
			# logger.debug(f'getting pose protein system {self}')
			self.protein_system = mp.parse(self.path, verbosity=False).protein_system
		return self._protein_system

	@protein_system.setter
	def protein_system(self, a):
		"""Sets the pose's protein molparse.System"""
		self._protein_system = a

	@property
	def complex_system(self):

		if self.has_complex_pdb_path:
			return mp.parse(self.path, verbosity=False)

		elif self.reference:

			# construct from .mol and reference

			system = self.reference.protein_system.copy()

			system.name = f'{self.target.name}_{self.reference.name}_{self.compound.name}'

			from molparse.rdkit import mol_to_AtomGroup
			ligand = mol_to_AtomGroup(self.mol)

			for atom in ligand.atoms:
				system.add_atom(atom)

			return system

		else:

			raise NotImplementedError

	@property
	def has_complex_pdb_path(self):
		return self.path.endswith('.pdb')

	@property
	def metadata(self):
		"""Returns the pose's metadata"""
		if self._metadata is None:
			self._metadata = self.db.get_metadata(table='pose', id=self.id)
		return self._metadata

	@property
	def fingerprint(self):
		"""Returns the pose's fingerprint"""
		return self._fingerprint

	@fingerprint.setter
	def fingerprint(self, fp):
		"""Set the pose's fingerprint"""

		# remove features that don't exist in this fingerprint?

		self._fingerprint = fp

		# store in the database
		fp = pickle.dumps(fp)
		# logger.debug('pickling fingerprint')
		self.db.update(table='pose', id=self.id, key=f'pose_fingerprint', value=fp)

	@property
	def tags(self):
		"""Returns the pose's tags"""
		if not self._tags:
			self._tags = self.get_tags()
		return self._tags

	@property
	def inspirations(self):
		"""Returns the pose's inspirations"""
		return self.get_inspirations()

	@property
	def features(self):
		"""Returns the pose's features"""
		return mp.rdkit.features_from_mol(self.mol)
	
	@property
	def dict(self):
		"""Serialised dictionary representing the pose"""
		return self.get_dict()

	@property
	def table(self):
		return self._table


	### METHODS

	def get_compound(self):
		return self.db.get_compound(id=self._compound_id)

	def get_tags(self):
		tags = self.db.select_where(query='tag_name', table='tag', key='pose', value=self.id, multiple=True, none='quiet')
		return TagSet(self, {t[0] for t in tags})
	
	def get_inspiration_ids(self):
		inspirations = self.db.select_where(query='inspiration_original', table='inspiration', key='derivative', value=self.id, multiple=True, none='quiet')
		if not inspirations:
			return None
		return set([v for v, in inspirations])

	def get_inspirations(self):	
		if not (inspirations := self.get_inspiration_ids()):
			return None

		from .pset import PoseSet
		return PoseSet(self.db, indices=inspirations)

	def get_dict(self, 
		mol: bool = False, 
		inspirations: bool | str = True, 
		reference: bool | str = True,
		metadata: bool = True,
		duplicate_name: str | bool = False,
	) -> dict:

		"""Returns a dictionary representing this Pose. Arguments:

		mol (bool): [True, False]
		inspirations (bool): [True, False, 'fragalysis']
		reference (bool): [True, False, 'name']

		"""

		serialisable_fields = ['id','name','smiles','path']

		data = {}
		for key in serialisable_fields:
			data[key] = getattr(self, key)

		if duplicate_name:
			assert isinstance(duplicate_name, str)
			data[duplicate_name] = data['name']

		if mol:
			data['mol'] = self.mol

		data['compound'] = self.compound.name
		data['target'] = self.target.name

		data['tags'] = self.tags
		
		if inspirations == 'fragalysis':
			data['inspirations'] = ','.join([p.name for p in self.inspirations])
		elif inspirations:
			data['inspirations'] = self.inspirations

		if reference == 'name':
			data['reference'] = self.reference.name
		elif reference:
			data['reference'] = self.reference

		if metadata and (metadict := self.metadata):
			for key in metadict:
				data[key] = metadict[key]

		return data

	def calculate_fingerprint(self):
		"""Calculate the pose's interaction fingerprint"""

		if self.path.endswith('.pdb'):
			
			import molparse as mp
			protein_system = self.protein_system
			if not self.protein_system:
				# logger.reading(self.path)
				protein_system = mp.parse(self.path, verbosity=False).protein_system

		elif self.path.endswith('.mol') and self.reference:

			logger.debug('fingerprint from .mol and reference pose')
			protein_system = self.reference.protein_system

		else:

			logger.debug('calculate_fingerprint()')
			raise NotImplementedError

		assert protein_system

		if not self.mol:
			return

		comp_features = self.features

		comp_features_by_family = {}
		for family in FEATURE_FAMILIES:
			comp_features_by_family[family] = [f for f in comp_features if f.family == family]

		# protein_features = self.target.features
		# if not protein_features:
		protein_features = self.target.calculate_features(protein_system)

		fingerprint = {}

		chains = protein_system.chain_names

		for prot_feature in protein_features:

			if prot_feature.chain_name not in chains:
				continue

			prot_family = prot_feature.family

			prot_residue = protein_system.get_chain(prot_feature.chain_name).residues[f'n{prot_feature.residue_number}']

			if not prot_residue:
				continue

			# if prot_feature.residue_number == 77:
			# 	logger.debug(repr(prot_feature))

			if prot_residue.name != prot_feature.residue_name:
				logger.warning(f'Feature {repr(prot_feature)}')
				continue

			prot_atoms = [prot_residue.get_atom(a).np_pos for a in prot_feature.atom_names.split(' ')]
			
			prot_coord = np.array(np.sum(prot_atoms,axis=0)/len(prot_atoms))

			complementary_family = COMPLEMENTARY_FEATURES[prot_family]

			complementary_comp_features = comp_features_by_family[complementary_family]

			cutoff = FEATURE_PAIR_CUTOFFS[f'{prot_family} {complementary_family}']

			valid_features = [f for f in complementary_comp_features if np.linalg.norm(f - prot_coord) <= cutoff]

			fingerprint[prot_feature.id] = len(valid_features)

		self.fingerprint = fingerprint

	def draw(self):
		"""Render this pose (and its inspirations)"""
		
		from molparse.rdkit import draw_mols
		
		mols = [self.mol]
		if self.inspirations:
			mols += [i.mol for i in self.inspirations]
		
		return draw_mols(mols)

	def grid(self):
		"""Draw a grid of this pose with its inspirations"""
		from molparse.rdkit import draw_grid
		from IPython.display import display

		mols = [self.compound.mol]
		labels = [self.plain_repr()]
		if self.inspirations:
			mols += [i.compound.mol for i in self.inspirations]
			labels += [i.plain_repr() for i in self.inspirations]

		display(draw_grid(mols, labels=labels))

	def summary(self, metadata:bool = True):
		"""Print a summary of this pose"""
		if self.name:
			logger.header(f'{str(self)}: {self.name}')
		else:
			logger.header(str(self))
		logger.var('smiles', self.smiles)
		logger.var('compound', repr(self.compound))
		logger.var('path', self.path)
		logger.var('target', repr(self.target))
		logger.var('reference', self.reference)
		logger.var('tags', self.tags)
		logger.var('inspirations', self.inspirations)
		if metadata:
			logger.var('metadata', str(self.metadata))

	def showcase(self):
		self.summary(metadata=False)
		self.grid()
		self.draw()
		from pprint import pprint
		logger.title('Metadata:')
		pprint(self.metadata)

	def plain_repr(self):
		"""Unformatted __repr__"""
		if self.name:
			return f'{self.compound}->{self}: {self.name}'
		else:
			return f'{self.compound}->{self}'

	### DUNDERS

	def __str__(self):
		return f'P{self.id}'

	def __repr__(self):
		return f'{mcol.bold}{mcol.underline}{self.plain_repr()}{mcol.unbold}{mcol.ununderline}'

	def __eq__(self, other):

		if isinstance(other, int):
			return self.id == other

		return self.id == other.id
