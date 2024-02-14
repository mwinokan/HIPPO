
import mcol

from rdkit import Chem

from .pose import Pose


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
	def poses(self):
		return self.get_poses()

	
	### METHODS

	def get_tags(self):
		tags = self.db.select_where(query='tag_name', table='tag', key='compound', value=self.id, multiple=True)
		return {t[0] for t in tags}

	def get_quotes(self, supplier=None):

		quote_ids = self.db.select_where(query='quote_id', table='quote', key='compound', value=self.id, multiple=True)

		if quote_ids:
			quotes = [self.db.get_quote(id=q[0]) for q in quote_ids]
		else:
			return []

		if supplier:
			quotes = [q for q in quotes if q['supplier'] == supplier]

		return quotes

	def get_poses(self, 
		target: str = None
	) -> list[Pose]:

		pose_ids = self.db.select_where(query='pose_id', table='pose', key='compound', value=self.id, multiple=True)

		poses = [self.db.get_pose(id=q[0]) for q in pose_ids]

		if target:
			poses = [q for q in poses if q['target'] == target]

		return poses


	### DUNDERS

	def __str__(self):
		return f'C{self.id}'

	def __repr__(self):
		# return f'Compound(#{self.id}, "{self.name}", {self.smiles})'
		# return f'[C{self.id} "{self.name}"]'
		return f'{mcol.bold}{mcol.underline}{self} "{self.name}"{mcol.unbold}{mcol.ununderline}'
		# return f'{mcol.bold}{mcol.underline}C{self.id}{mcol.unbold}{mcol.ununderline}'

