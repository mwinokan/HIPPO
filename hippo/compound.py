
import mcol

from rdkit import Chem

from .pose import Pose
from .tags import TagSubset

import logging
logger = logging.getLogger('HIPPO')


class Compound:

	def __init__(self,
			db,
			id: int,
			name: str,
			smiles: str,
			base: int,
			mol: Chem.Mol | bytes | None = None,
			metadata: dict | None = None,
	):
		
		self._id = id
		self._name = name
		self._smiles = smiles
		self._base = base
		self._table = 'compound'
		self._tags = None

		self._metadata = metadata
		
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
		if self._mol is None:
			mol, = self.db.select_where(query='mol_to_binary_mol(compound_mol)', table='compound', key='id', value=self.id, multiple=False)
			self._mol = Chem.Mol(mol)
		return self._mol

	@property
	def metadata(self):
		if self._metadata is None:
			self._metadata = self.db.get_metadata(table='compound', id=self.id)
		return self._metadata

	@property
	def db(self):
		return self._db

	@property
	def tags(self):
		if not self._tags:
			self._tags = self.get_tags()
		return self._tags

	@property
	def poses(self):
		return self.get_poses()

	@property
	def num_poses(self):
		return self.db.count_where(table='pose', key='compound', value=self.id)

	@property
	def base(self):
		if isinstance(self._base, int):
			self._base = self.db.get_compound(id=self._base)
		return self._base

	@base.setter
	def base(self, b):
		self.set_base(b)

	@property
	def table(self):
		return self._table

	@property
	def reactions(self):
		return self.get_reactions(none=False)
	
	@property
	def dict(self):

		serialisable_fields = ['id','name','smiles']

		data = {}
		for key in serialisable_fields:
			data[key] = getattr(self, key)

		if base := self.base:
			data['base'] = base.name

		if metadata := self.metadata:
			for key in metadata:
				data[key] = metadata[key]

		return data
	
	### METHODS

	def get_tags(self) -> set:
		tags = self.db.select_where(query='tag_name', table='tag', key='compound', value=self.id, multiple=True, none='quiet')
		return TagSubset(self, {t[0] for t in tags}, commit=False)
		# if tags:
		# 	return TagSubset(self, {t[0] for t in tags})
		# else:
		# 	return None

	def get_quotes(self, min_amount=None, supplier=None, max_lead_time=None, none='error', pick_cheapest=False, df=False) -> list[dict]:

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

		if df:
			from pandas import DataFrame
			return DataFrame([q.asdict() for q in quotes]).drop(columns='compound')
		
		return quotes

	def get_reactions(self, none='error') -> list:

		reaction_ids = self.db.select_where(query='reaction_id', table='reaction', key='product', value=self.id, multiple=True, none=none)

		if reaction_ids:
			reactions = [self.db.get_reaction(id=q[0]) for q in reaction_ids]
		else:
			return []

		if len(reactions) > 1:
			reactions = self.db.prune_reactions(compound=self, reactions=reactions)

		return reactions

	def get_poses(self, 
		target: str = None
	) -> list[Pose]:

		pose_ids = self.db.select_where(query='pose_id', table='pose', key='compound', value=self.id, multiple=True)

		poses = [self.db.get_pose(id=q[0]) for q in pose_ids]

		if target:
			poses = [q for q in poses if q['target'] == target]

		return poses

	def set_base(self, base, commit=True):
		self._base = base
		self.db.update(table='compound', id=self.id, key='compound_base', value=base, commit=commit)

	def as_ingredient(self, amount, max_lead_time=None, supplier=None):
		
		quote = self.get_quotes(
			pick_cheapest=True, 
			min_amount=amount, 
			max_lead_time=max_lead_time, 
			supplier=supplier,
			none='quiet',
		)
		
		return Ingredient(self, amount, quote, max_lead_time, supplier)

	def draw(self):

		if self.base:
			from molparse.rdkit import draw_mcs
			return draw_mcs({self.base.smiles:f'{self.base} (base)', self.smiles:str(self)})
		else:
			display(self.mol)

	### DUNDERS

	def __str__(self):
		return f'C{self.id}'

	def __repr__(self):
		return f'{mcol.bold}{mcol.underline}{self} "{self.name}"{mcol.unbold}{mcol.ununderline}'

	def __eq__(self, other):
		return self.id == other.id

class Ingredient(Compound):

	"""An ingredient is a Compound with a fixed quanitity and an attached quote"""

	def __init__(self, inherit, amount, quote, max_lead_time=None, supplier=None):
		self._id = inherit.id
		self._name = inherit.name
		self._smiles = inherit.smiles
		self._base = inherit.base			
		self._mol = inherit.mol
		self._db = inherit.db
		
		self._max_lead_time = max_lead_time
		self._supplier = supplier
		self._amount = amount
		self._quote = quote

	@property
	def amount(self):
		return self._amount

	@amount.setter
	def amount(self, a):

		quote = self.get_quotes(
			pick_cheapest=True, 
			min_amount=a, 
			max_lead_time=self._max_lead_time, 
			supplier=self._supplier,
			none='quiet',
		)
		
		self._amount = a

	@property
	def quote(self):
		return self._quote

	def __repr__(self):
		return f'{self.amount:.2f}mg of {mcol.bold}{mcol.underline}{self} "{self.name}"{mcol.unbold}{mcol.ununderline}'

	def __eq__(self, other):

		if self.id != other.id:
			return False

		return self.amount == other.amount
