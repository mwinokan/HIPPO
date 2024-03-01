
import mcol
import molparse as mp

from rdkit import Chem

from .tags import TagSet

import logging
logger = logging.getLogger('HIPPO')

class Pose:

	def __init__(self, 
		db,
		id: int,
		name: str,
		longname: str,
		smiles: str,
		reference: int, # another pose
		path: str,
		compound: int,
		target: str,
		mol: Chem.Mol | bytes | None,
		fingerprint: str,
		metadata: dict | None = None,
	):

		self._db = db
		self._id = id
		self._name = name
		self._longname = longname
		self._smiles = smiles
		self._compound_id = compound
		self._target = target
		self._path = path
		self._fingerprint = fingerprint
		self._metadata = metadata
		self._tags = None
		self._table = 'pose'
		self._reference = reference

		if isinstance(mol, bytes):
			self._mol = Chem.Mol(mol)
		else: 
			self._mol = mol
		
	### FACTORIES

	### PROPERTIES

	@property
	def db(self):
		return self._db

	@property
	def id(self):
		return self._id

	@property
	def name(self):
		return self._name

	@property
	def longname(self):
		return self._longname

	@property
	def smiles(self):
		return self._smiles

	@property
	def target(self):
		if isinstance(self._target, int):
			self._target = self.db.get_target(id=self._target)
		return self._target

	@property
	def compound_id(self):
		return self._compound_id

	@property
	def compound(self):
		return self.get_compound()

	@property
	def path(self):
		return self._path

	@property
	def reference(self):
		if isinstance(self._reference, int):
			self._reference = self.db.get_pose(id=self._reference)
		return self._reference

	@property
	def mol(self):
		if not self._mol and self.path:
			logger.reading(self.path)
			sys = mp.parse(self.path, verbosity=False)
			lig_residues = sys['rLIG']
			if len(lig_residues) > 1:
				logger.warning('Multiple ligands in PDB')
			lig_res = lig_residues[0]
			self.mol = lig_res.rdkit_mol
		return self._mol

	@mol.setter
	def mol(self, m):

		self._mol = m

		self.db.update(table='pose', id=self.id, key='pose_mol', value=m.ToBinary())

	@property
	def metadata(self):
		if self._metadata is None:
			self._metadata = self.db.get_metadata(table='pose', id=self.id)
		return self._metadata

	@property
	def fingerprint(self):
		return self._fingerprint

	@property
	def tags(self):
		if not self._tags:
			self._tags = self.get_tags()
		return self._tags

	@property
	def inspirations(self):
		return self.get_inspirations()

	@property
	def table(self):
		return self._table
	

	### METHODS

	def get_compound(self):
		return self.db.get_compound(id=self._compound_id)

	def get_tags(self):
		tags = self.db.select_where(query='tag_name', table='tag', key='pose', value=self.id, multiple=True)
		return TagSet(self, {t[0] for t in tags})
	
	def get_inspirations(self):
		inspirations = self.db.select_where(query='inspiration_original', table='inspiration', key='derivative', value=self.id, multiple=True, none='quiet')
		
		if inspirations:
			inspirations = [self.db.get_pose(id=id[0]) for id in inspirations]

		return inspirations

	### DUNDERS

	def __str__(self):
		return f'P{self.id}'

	def __repr__(self):
		return f'{mcol.bold}{mcol.underline}{self.compound}->{self} "{self.longname}"{mcol.unbold}{mcol.ununderline}'
		