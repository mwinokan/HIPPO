
# from .tools import df_row_to_dict

from .compound import Compound
from .db import Database

from .recipe import Recipe

import mcol

import os

import logging
logger = logging.getLogger('HIPPO')

class CompoundSet:
	"""Object representing the 'compound' table in the :class:`.Database`."""

	def __init__(self, 
		db: Database, 
		table: str = 'compound',
	) -> None:
		
		self._db = db
		self._table = table

	### FACTORIES

	### PROPERTIES

	@property
	def db(self) -> Database:
		"""Returns the associated :class:`.database`"""
		return self._db
	
	@property
	def table(self) -> str:
		return self._table

	@property
	def names(self) -> list:
		"""Returns the names of child compounds"""
		result = self.db.select(table=self.table, query='compound_name', multiple=True)
		return [q for q, in result]

	@property
	def ids(self):
		"""Returns the IDs of child compounds"""
		result = self.db.select(table=self.table, query='compound_id', multiple=True)
		return [q for q, in result]
	
	@property
	def tags(self):
		"""Returns the set of unique tags present in this compound set"""
		values = self.db.select_where(table='tag', query='DISTINCT tag_name', key='tag_compound IS NOT NULL', multiple=True)
		return set(v for v, in values)

	### METHODS

	def get_by_tag(self,tag):
		"""Get all child compounds with a certain tag"""
		values = self.db.select_where(query='tag_compound', table='tag', key='name', value=tag, multiple=True)
		ids = [v for v, in values if v]
		return self[ids]

	def get_by_metadata(self, key: str, value: str | None = None):
		"""Get all child compounds with by their metadata. If no value is passed, then simply containing the key in the metadata dictionary is sufficient"""
		results = self.db.select(query='compound_id, compound_metadata', table='compound', multiple=True)
		if value is None:
			ids = [i for i,d in results if d and f'"{key}":' in d]
		else:
			if isinstance(value, str):
				value = f'"{value}"'
			ids = [i for i,d in results if d and f'"{key}": {value}' in d]
		return self[ids]		

	def summary(self):
		"""Print a summary of this compound set"""
		logger.header('CompoundSet()')
		logger.var('#compounds', len(self))
		# logger.var('#poses', self.num_poses)
		logger.var('tags', self.tags)

	### DUNDERS

	def __getitem__(self, key) -> Compound:
		
		match key:

			case int():

				if key == 0:
					return self.__getitem__(key=1)

				if key < 0:
					key = len(self) + 1 + key
					return self.__getitem__(key=key)

				else:
					return self.db.get_compound(table=self.table, id=key)

			case str():
				return self.db.get_compound(table=self.table, name=key)

			case list():
				return CompoundSubset(self.db, key)

			case tuple():
				return CompoundSubset(self.db, key)

			case slice():

				start = key.start or 1
				stop = key.stop or len(self)
				step = key.step or 1

				indices = [i for i in range(start, stop, step)]

				return CompoundSubset(self.db, indices)

			case _:
				logger.error(f'Unsupported type for CompoundSet.__getitem__(): {key=} {type(key)}')

		return None

	def __repr__(self) -> str:
		return f'{mcol.bold}{mcol.underline}set(C x {len(self)}){mcol.unbold}{mcol.ununderline}'

	def __len__(self) -> int:
		return self.db.count(self.table)

	def __iter__(self):
		return iter(self[i+1] for i in range(len(self)))


class CompoundSubset(CompoundSet):
	"""Object representing a subset of the 'compound' table in the :class:`.Database`."""

	def __init__(self,
		db: Database,
		indices: list = None,
		*,
		table: str = 'compound',
	):

		self._db = db
		self._table = table

		self._indices = indices or []

		self._indices = list(set(self.indices))

	### PROPERTIES

	@property
	def indices(self):
		"""Returns the ids of compounds in this set"""
		return self._indices

	@property
	def ids(self):
		"""Returns the ids of compounds in this set"""
		return self.indices

	@property
	def names(self):
		"""Returns the names of compounds in this set"""
		return [self.db.select_where(table=self.table, query='compound_name', key='id', value=i, multiple=False)[0] for i in self.indices]

	@property
	def tags(self):
		"""Returns the set of unique tags present in this compound set"""
		values = self.db.select_where(table='tag', query='DISTINCT tag_name', key=f'tag_compound in {tuple(self.ids)}', multiple=True)
		return set(v for v, in values)

	@property
	def num_poses(self):
		"""Count the poses associated to this set of compounds"""
		from .pset import PoseSubset
		return self.db.count_where(table='pose', key=f'pose_compound in {tuple(self.ids)}')

	@property
	def poses(self):
		"""Get the poses associated to this set of compounds"""
		from .pset import PoseSubset
		ids = self.db.select_where(query='pose_id', table='pose', key=f'pose_compound in {tuple(self.ids)}', multiple=True)
		ids = [v for v, in ids]
		return PoseSubset(self.db, ids)

	### METHODS

	def get_by_tag(self,tag):
		"""Get all child compounds with a certain tag"""
		values = self.db.select_where(query='tag_compound', table='tag', key='name', value=tag, multiple=True)
		ids = [v for v, in values if v and v in self.ids]
		return CompoundSubset(self.db, ids)

	def get_by_metadata(self, key: str, value: str | None = None):
		"""Get all child compounds with by their metadata. If no value is passed, then simply containing the key in the metadata dictionary is sufficient"""
		results = self.db.select(query='compound_id, compound_metadata', table='compound', multiple=True)
		if value is None:
			ids = [i for i,d in results if d and f'"{key}":' in d and i in self.ids]
		else:
			if isinstance(value, str):
				value = f'"{value}"'
			ids = [i for i,d in results if d and f'"{key}": {value}' in d and i in self.ids]
		return CompoundSubset(self.db, ids)		

	def draw(self):
		"""Draw a grid of all contained molecules"""
		from molparse.rdkit import draw_grid

		data = [(str(c), c.mol) for c in self]

		mols = [d[1] for d in data]
		labels = [d[0] for d in data]

		return draw_grid(mols, labels=labels)

	def get_recipe(self, 
		amount: float = 1, 
		simplify_network: bool = True,
	):
		"""Generate the :class:`.Recipe` to make these compounds."""

		n_comps = len(self)

		assert simplify_network

		if not hasattr(amount, '__iter__'):
			amount = [amount] * n_comps

		products = []
		reactants = []
		reactions = []
		intermediates = []

		for comp, a in zip(self, amount):

			reax = comp.reactions

			assert len(reax) == 1

			recipe = reax[0].get_recipe(a)

			products.append(comp.as_ingredient(amount=a))

			for reactant in recipe.reactants:
				matches = [r for r in reactants if r.id == reactant.id]
				if matches:
					matches[0].amount += reactant.amount
				else:
					reactants.append(reactant) 

			for reaction in recipe.reactions:
				matches = [r for r in reactions if r.id == reaction.id]
				if not matches:
					reactions.append(reaction) 

			for intermediate in recipe.intermediates:
				matches = [r for r in intermediates if r.id == intermediate.id]
				if not matches:
					intermediates.append(intermediate) 

		return Recipe(products=products, reactants=reactants, reactions=reactions, intermediates=intermediates)

	def summary(self):
		"""Print a summary of this compound set"""
		logger.header('CompoundSubset()')
		logger.var('#compounds', len(self))
		logger.var('#poses', self.num_poses)
		logger.var('tags', self.tags)

	### DUNDERS

	def __repr__(self) -> str:
		return f'{mcol.bold}{mcol.underline}subset(C x {len(self)}){mcol.unbold}{mcol.ununderline}'

	def __len__(self) -> int:
		return len(self.indices)

	def __iter__(self):
		return iter(self.db.get_compound(table=self.table, id=i) for i in self.indices)

	def __getitem__(self, key) -> Compound:
		try:
			index = self.indices[key]
		except IndexError:
			logger.exception(f'list index out of range: {key=} for {self}')
			raise
		return self.db.get_pose(table=self.table, id=index)
