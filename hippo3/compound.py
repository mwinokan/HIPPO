
from rdkit import Chem


class Compound:

	def __init__(self,
			db,
			id: int,
			name: str,
			smiles: str,
			base: int,
			mol: Chem.Mol,
		):
		
		self._id = id
		self._name = name
		self._smiles = smiles
		self._base = base
		
		if isinstance(mol, bytes):
			mol = Chem.Mol(mol)
			
		self._mol = mol

		self._db = db
		
	### FACTORIES

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
	def db(self):
		return self._db

	@property
	def tags(self):
		return self.get_tags()

	@property
	def inspirations(self):
		return self.get_inspirations()
	
	### METHODS

	def get_tags(self):
		tags = self.db.select_where(query='tag_name', table='tag', key='compound', value=self.id, multiple=True)
		return {t[0] for t in tags}

	def get_inspirations(self):
		inspirations = self.db.select_where(query='inspiration_pose', table='inspiration', key='compound', value=self.id, multiple=True, none='quiet')
		
		if inspirations:
			inspirations = [self.db.get_pose(id=id[0]) for id in inspirations]

		return inspirations

	def get_quotes(self, supplier=None):

		quote_ids = self.db.select_where(query='quote_id', table='quote', key='compound', value=self.id, multiple=True)

		quotes = [self.db.get_quote(id=q[0]) for q in quote_ids]

		if supplier:
			quotes = [q for q in quotes if q['supplier'] == supplier]

		return quotes

	### DUNDERS

	def __repr__(self):
		return f'Compound({self.name}, {self.smiles})'
