
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

	@property
	def base(self):
		return self._base

	
	### METHODS

	def get_tags(self) -> set:
		tags = self.db.select_where(query='tag_name', table='tag', key='compound', value=self.id, multiple=True)
		return {t[0] for t in tags}

	def get_quotes(self, min_amount=None, supplier=None, max_lead_time=None, none='error', pick_cheapest=False) -> list[dict]:

		quote_ids = self.db.select_where(query='quote_id', table='quote', key='compound', value=self.id, multiple=True, none=none)

		if quote_ids:
			quotes = [self.db.get_quote(id=q[0]) for q in quote_ids]
		else:
			return []

		if supplier:
			quotes = [q for q in quotes if q.supplier == supplier]

		if min_amount:
			quotes = [q for q in quotes if q.amount >= min_amount]

		if max_lead_time:
			quotes = [q for q in quotes if q.lead_time <= max_lead_time]

		if pick_cheapest:
			return sorted(quotes, key=lambda x: x.price)[0]
		
		return quotes

	def get_reactions(self, none='error') -> list:

		reaction_ids = self.db.select_where(query='reaction_id', table='reaction', key='product', value=self.id, multiple=True, none=none)

		if reaction_ids:
			reactions = [self.db.get_reaction(id=q[0]) for q in reaction_ids]
		else:
			return []

		return reactions

	def get_poses(self, 
		target: str = None
	) -> list[Pose]:

		pose_ids = self.db.select_where(query='pose_id', table='pose', key='compound', value=self.id, multiple=True)

		poses = [self.db.get_pose(id=q[0]) for q in pose_ids]

		if target:
			poses = [q for q in poses if q['target'] == target]

		return poses

	def as_ingredient(self, amount, max_lead_time=None, supplier=None):
		
		quote = self.get_quotes(
			pick_cheapest=True, 
			min_amount=amount, 
			max_lead_time=max_lead_time, 
			supplier=supplier
		)
		
		return Ingredient(self, amount, quote)

	### DUNDERS

	def __str__(self):
		return f'C{self.id}'

	def __repr__(self):
		return f'{mcol.bold}{mcol.underline}{self} "{self.name}"{mcol.unbold}{mcol.ununderline}'

class Ingredient(Compound):

	"""An ingredient is a Compound with a fixed quanitity and an attached quote"""

	def __init__(self, inherit, amount, quote):
		self._id = inherit.id
		self._name = inherit.name
		self._smiles = inherit.smiles
		self._base = inherit.base			
		self._mol = inherit.mol
		self._db = inherit.db

		self._required_amount = amount
		self._quote = quote

	@property
	def required_amount(self):
		return self._required_amount

	@property
	def quote(self):
		return self._quote

	def __repr__(self):
		return f'{self.required_amount:.2f}mg of {mcol.bold}{mcol.underline}{self} "{self.name}"{mcol.unbold}{mcol.ununderline}'