
# from .tools import df_row_to_dict

from .compound import Compound
from .db import Database

from .recipe import Recipe

import mcol

import os

import logging
logger = logging.getLogger('HIPPO')

class CompoundSet:

	# class to access entries in database tables containing compounds

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
		return self._db
	
	@property
	def table(self) -> str:
		return self._table

	@property
	def names(self):
		result = self.db.select(table=self.table, query='compound_name', multiple=True)
		return [q for q, in result]

	@property
	def ids(self):
		result = self.db.select(table=self.table, query='compound_id', multiple=True)
		return [q for q, in result]
	
	### METHODS

	def get_by_tag(self,tag):
		values = self.db.select_where(query='tag_compound', table='tag', key='name', value=tag, multiple=True)
		ids = [v for v, in values if v]
		return self[ids]

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
				return CompoundSubset(self.db, self.table, key)

			case tuple():
				return CompoundSubset(self.db, self.table, key)

			case slice():

				start = key.start or 1
				stop = key.stop or len(self)
				step = key.step or 1

				indices = [i for i in range(start, stop, step)]

				return CompoundSubset(self.db, self.table, indices)

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

	def __init__(self,
		db: Database,
		table: str = 'compound',
		indices: list = None,
	):

		self._db = db
		self._table = table

		self._indices = indices or []

		self._indices = list(set(self.indices))

	### PROPERTIES

	@property
	def indices(self):
		return self._indices

	@property
	def ids(self):
		return self.indices

	@property
	def names(self):
		return [self.db.select_where(table=self.table, query='compound_name', key='id', value=i, multiple=False)[0] for i in self.indices]

	### METHODS

	def get_by_tag(self,tag):
		values = self.db.select_where(query='tag_compound', table='tag', key='name', value=tag, multiple=True)
		ids = [v for v, in values if v]
		ids = [v for v in ids if v in self.indices]
		return CompoundSubset(self.db, self.table, ids)

	def draw(self):
		from molparse.rdkit import draw_grid

		data = [(str(c), c.mol) for c in self]

		mols = [d[1] for d in data]
		labels = [d[0] for d in data]

		return draw_grid(mols, labels=labels)

	def get_recipe(self, amount=1):

		n_comps = len(self)

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
