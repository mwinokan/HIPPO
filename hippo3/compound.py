
from rdkit import Chem

class Compound:

	def __init__(self,
			db,
			id,
			name,
			smiles,
			base,
			mol,
		):
		
		self._id = id
		self._name = name
		self._smiles = smiles
		self._base = base
		self._mol = mol
		self._db = db
		
	### FACTORIES

	@classmethod
	def from_db_entry(cls, 
		db,
		id,
		name,
		smiles,
		base,
		binary_mol,
	):

		self = cls.__new__(cls)

		mol = Chem.Mol(binary_mol)

		self.__init__(db, id, name, smiles, base, mol)
		
		return self

	### PROPERTIES

	@property
	def id(self):
		return self._id
	
	@property
	def name(self):
		return self._name
	
	@property
	def smiles(self):
		return self._smiles
	
	@property
	def mol(self):
		return self._mol

	@property
	def tags(self):
		return self.get_tags()

	@property
	def db(self):
		return self._db
	
	### METHODS

	def get_tags(self):
	
		tags = self.db.select_where(query='tag_name', table='tag', key='compound', value=self.id, multiple=True)

		return set({t[0] for t in tags})

	### DUNDERS

	def __repr__(self):
		return f'Compound({self.name}, {self.smiles})'
