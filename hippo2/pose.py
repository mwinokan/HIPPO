
import mout

from .fingerprint import Fingerprint

class Pose:

	def __init__(self, name, compound, pdb_path, site_index, chain=None, tags=None):

		# arguments
		self._name = name
		self._compound = compound
		self._pdb_path = pdb_path
		self.site_index = site_index
		self._chain = chain

		self.tags = list(tags or [])

		# blanks
		self._fingerprint = None
		self._pdb_entry = None
		self._mol = None
		
	### FACTORIES

	@classmethod
	def from_bound_pdb(cls, compound, pdb_path, metadata, tags=None, site_index=None, chain=None):

		self = cls.__new__(cls)

		tags = list(tags or [])

		if metadata['site_name']:
			tags.append(metadata['site_name'])
		
		if metadata['pdb_entry']:
			self._pdb_entry = metadata['pdb_entry']
			self.tags.append('InPDB')

		pose_name = metadata['crystal_name']

		if site_index is None:
			site_index = pose_name[-2]

		if chain is None:
			chain = pose_name[-1]

		self.__init__(f'{site_index}{chain}', compound, pdb_path, site_index, chain, tags)
		
		return self

	### PROPERTIES

	@property
	def compound(self):
		return self._compound

	@property
	def name(self):
		return self._name

	@property
	def longname(self):
		return f'{self.compound.name}_{self.name}'
	
	@property
	def site_index(self):
		return self._site_index

	@site_index.setter
	def site_index(self, i):
		self._site_index = int(i)
	
	@property
	def chain(self):
		return self._chain

	@property
	def pdb_path(self):
		return self._pdb_path
	
	@property
	def fingerprint(self):

		if self._fingerprint is None:
			self._fingerprint = self.calculate_fingerprint()

		return self._fingerprint

	@property
	def bound_system(self):
		import molparse as mp
		sys = mp.parsePDB(self.pdb_path, dry=True, verbosity=0)
		return sys

	@property
	def mol(self):
		
		if self._mol is None:
		
			lig_residues = self.bound_system['rLIG']
			lig_residues = [l for l in lig_residues if l.chain == self.chain]
			
			split_lig_residues = []
			for lig in lig_residues:
				split_lig_residues += lig.split_by_site()
			
			ligand_group = split_lig_residues[self.site_index]
			
			from molparse.rdkit import mol_from_pdb_block
			
			self._mol = mol_from_pdb_block(ligand_group.pdb_block)
			
		return self._mol

	### METHODS
	
	def calculate_fingerprint(self):

		# need to load the bound PDB
		sys = self.bound_system

		fingerprint = Fingerprint(self, self.mol, sys.protein_system)
		
		# fingerprint = None
		return fingerprint

	### DUNDERS

	def __repr__(self):
		return f'Pose("{self.name}", compound={self.compound.name}, tags={self.tags})'

