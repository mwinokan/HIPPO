
import mcol

from rdkit import Chem

from .pose import Pose
# from .pset import PoseSubset
from .tags import TagSubset
# from .rset import ReactionSet

import logging
logger = logging.getLogger('HIPPO')


class Compound:

	"""A :class:`.Compound` represents a ligand/small molecule with stereochemistry removed and no atomic coordinates. I.e. it represents the chemical structure. It's name is always an InChiKey. If a compound is an elaboration it can have a :meth:`.Compound.base` property which is another :class:`.Compound`. :class:`.Compound` objects are target-agnostic and can be linked to any number of catalogue entries (:class:`.Quote`) or synthetic pathways (:class:`.Reaction`).

	:class:`.Compound` objects should not be created directly. Instead use :meth:`.HIPPO.register_compound` or :meth:`.HIPPO.compounds`
	"""

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
		self._tags = None
		self._table = 'compound'

		self._metadata = metadata
		
		if isinstance(mol, bytes):
			mol = Chem.Mol(mol)
			
		self._mol = mol

		self._db = db
		
	### FACTORIES

	### PROPERTIES

	@property
	def id(self) -> int:
		"""Returns the compound's database ID"""
		return self._id
	
	@property
	def name(self) -> str:
		"""Returns the compound's name (InChiKey)"""
		return self._name
	
	@property
	def smiles(self) -> str:
		"""Returns the compound's (flattened) smiles"""
		return self._smiles
	
	@property
	def mol(self) -> Chem.Mol:
		"""Returns the compound's RDKit Molecule"""
		if self._mol is None:
			mol, = self.db.select_where(query='mol_to_binary_mol(compound_mol)', table='compound', key='id', value=self.id, multiple=False)
			self._mol = Chem.Mol(mol)
		return self._mol

	@property
	def metadata(self) -> dict:
		"""Returns the compound's metadata dict"""
		if self._metadata is None:
			self._metadata = self.db.get_metadata(table='compound', id=self.id)
		return self._metadata

	@property
	def db(self):
		"""Returns a pointer to the parent database"""
		return self._db

	@property
	def tags(self) -> TagSubset:
		"""Returns the compound's tags"""
		if not self._tags:
			self._tags = self.get_tags()
		return self._tags

	@property
	def poses(self):
		"""Returns the compound's poses"""
		return self.get_poses()

	@property
	def num_poses(self) -> int:
		"""Returns the number of associated poses"""
		return self.db.count_where(table='pose', key='compound', value=self.id)

	@property
	def num_reactions(self) -> int:
		"""Returns the number of associated reactions (product)"""
		return self.db.count_where(table='reaction', key='product', value=self.id)

	@property
	def num_reactant(self) -> int:
		"""Returns the number of associated reactions (reactant)"""
		return self.db.count_where(table='reactant', key='compound', value=self.id)

	@property
	def base(self):
		"""Returns the base compound for this elaboration"""
		if isinstance(self._base, int):
			self._base = self.db.get_compound(id=self._base)
		return self._base

	@base.setter
	def base(self, b):
		"""Set the base compound for this elaboration"""
		self.set_base(b)

	@property
	def reactions(self):
		"""Returns the reactions resulting in this compound"""
		return self.get_reactions(none=False, prune_duplicate=False)

	@property
	def dict(self) -> dict:
		"""Returns a dictionary of this compound"""

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

	@property
	def table(self):
		return self._table
	
	
	### METHODS

	def get_tags(self) -> set:
		tags = self.db.select_where(query='tag_name', table='tag', key='compound', value=self.id, multiple=True, none='quiet')
		return TagSubset(self, {t[0] for t in tags}, commit=False)

	def get_quotes(self, min_amount=None, supplier=None, max_lead_time=None, none='error', pick_cheapest=False, df=False) -> list[dict]:
		"""Get all quotes associated to this compound"""

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

	def get_reactions(self, none='error', as_reactant=False, prune_duplicate=True) -> list:
		"""Get the associated reactions as product, unless as_reactant is True."""

		if as_reactant:
			reaction_ids = self.db.select_where(query='reactant_reaction', table='reactant', key='compound', value=self.id, multiple=True, none=none)
		else:
			reaction_ids = self.db.select_where(query='reaction_id', table='reaction', key='product', value=self.id, multiple=True, none=none)

		if reaction_ids:
			reactions = [self.db.get_reaction(id=q[0]) for q in reaction_ids]
		else:
			return []

		if not as_reactant and prune_duplicate and len(reactions) > 1:
			reactions = self.db.prune_reactions(compound=self, reactions=reactions)

		return reactions

	def get_poses(self, 
		target: str = None
	) -> list[Pose]:

		pose_ids = self.db.select_where(query='pose_id', table='pose', key='compound', value=self.id, multiple=True)

		# poses = [self.db.get_pose(id=q[0]) for q in pose_ids]

		from .pset import PoseSubset

		return PoseSubset(self.db, [q[0] for q in pose_ids])

	def set_base(self, base, commit=True):
		self._base = base
		self.db.update(table='compound', id=self.id, key='compound_base', value=base, commit=commit)

	def as_ingredient(self, amount, max_lead_time=None, supplier=None):
		"""Convert this compound into an Ingredient object"""
		
		quote = self.get_quotes(
			pick_cheapest=True, 
			min_amount=amount, 
			max_lead_time=max_lead_time, 
			supplier=supplier,
			none='quiet',
		)
		
		return Ingredient(self, amount, quote, max_lead_time, supplier)

	def draw(self):
		"""Display this compound (and its base if it has one)"""

		if self.base:
			from molparse.rdkit import draw_mcs
			return draw_mcs({self.base.smiles:f'{self.base} (base)', self.smiles:str(self)})
		else:
			display(self.mol)

	def summary(self):
		"""Print a summary of this compound"""
		logger.header(repr(self))
		logger.var('smiles', self.smiles)
		logger.var('base', self.base)
		logger.var('#poses', self.num_poses)
		logger.var('#reactions (product)', self.num_reactions)
		logger.var('#reactions (reactant)', self.num_reactant)
		logger.var('tags', self.tags)
		logger.var('metadata', str(self.metadata))

	### DUNDERS

	def __str__(self):
		return f'C{self.id}'

	def __repr__(self):
		return f'{mcol.bold}{mcol.underline}{self} "{self.name}"{mcol.unbold}{mcol.ununderline}'

	def __eq__(self, other):
		return self.id == other.id

class Ingredient(Compound):

	"""An ingredient is a :class:`.Compound` with a fixed quanitity and an attached quote.

	Create one from a :meth:`.Compound.as_ingredient`"""

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
	def amount(self) -> float:
		"""Returns the amount"""
		return self._amount

	@amount.setter
	def amount(self, a):
		"""Set the amount"""

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
		"""Returns the associated :class:`Quote`"""
		return self._quote

	def __repr__(self):
		return f'{self.amount:.2f}mg of {mcol.bold}{mcol.underline}{self} "{self.name}"{mcol.unbold}{mcol.ununderline}'

	def __eq__(self, other):

		if self.id != other.id:
			return False

		return self.amount == other.amount
