
# from .tools import df_row_to_dict

from .compound import Compound, Ingredient
from .db import Database

from .recipe import Recipe

import mcol

import os

import logging
logger = logging.getLogger('HIPPO')

class CompoundTable:
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

	def get_by_base(self,base):

		if not isinstance(base, int):
			assert base._table == 'compound'
			base = base.id

		values = self.db.select_where(query='compound_id', table='compound', key='base', value=base, multiple=True)
		ids = [v for v, in values if v]
		return self[ids]
	
	def summary(self):
		"""Print a summary of this compound set"""
		logger.header('CompoundTable()')
		logger.var('#compounds', len(self))
		# logger.var('#poses', self.num_poses)
		logger.var('tags', self.tags)

	def draw(self):
		return self[self.ids].draw()

	def interactive(self):
		return self[self.ids].interactive()

	def get_df(self, **kwargs):

		from pandas import DataFrame

		data = []

		for pose in self:
			d = pose.get_dict(**kwargs)
			data.append(d)

		return DataFrame(data)

	### DUNDERS

	def __call__(self, *, tag=None, base=None):
		if tag:
			return self.get_by_tag(tag)
		elif base:
			return self.get_by_base(base)
		else:
			raise NotImplementedError

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
				comp = self.db.get_compound(table=self.table, inchikey=key)
				if not comp:
					comp = self.db.get_compound(table=self.table, alias=key)
				return comp

			case key if isinstance(key, list) or isinstance(key, tuple):

				indices = []
				for i in key:
					if isinstance(i,int):
						index = i
					elif isinstance(i,str):
						index = self.db.get_compound_id(name=i)
					else:
						raise NotImplementedError

					assert index
					indices.append(index)

				return CompoundSet(self.db, indices)

			case slice():

				start = key.start or 1
				stop = key.stop or len(self)
				step = key.step or 1

				indices = [i for i in range(start, stop, step)]

				return CompoundSet(self.db, indices)

			case _:
				logger.error(f'Unsupported type for CompoundTable.__getitem__(): {key=} {type(key)}')

		return None

	def __repr__(self) -> str:
		return f'{mcol.bold}{mcol.underline}set(C x {len(self)}){mcol.unbold}{mcol.ununderline}'

	def __len__(self) -> int:
		return self.db.count(self.table)

	def __iter__(self):
		return iter(self[i+1] for i in range(len(self)))


class CompoundSet(CompoundTable):
	"""Object representing a subset of the 'compound' table in the :class:`.Database`."""

	def __init__(self,
		db: Database,
		indices: list = None,
		*,
		table: str = 'compound',
	):

		self._db = db
		self._table = table

		indices = indices or []

		if not isinstance(indices, list):
			indices = list(indices)

		assert all(isinstance(i, int) for i in indices)

		self._indices = sorted(list(set(indices)))

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
		"""Returns the aliases of compounds in this set"""
		return [self.db.select_where(table=self.table, query='compound_alias', key='id', value=i, multiple=False)[0] for i in self.indices]

	@property
	def inchikeys(self):
		"""Returns the inchikeys of compounds in this set"""
		return [self.db.select_where(table=self.table, query='compound_inchikey', key='id', value=i, multiple=False)[0] for i in self.indices]

	@property
	def tags(self):
		"""Returns the set of unique tags present in this compound set"""
		values = self.db.select_where(table='tag', query='DISTINCT tag_name', key=f'tag_compound in {tuple(self.ids)}', multiple=True)
		if not values:
			return set()
		return set(v for v, in values)

	@property
	def num_poses(self):
		"""Count the poses associated to this set of compounds"""
		from .pset import PoseSet
		return self.db.count_where(table='pose', key=f'pose_compound in {tuple(self.ids)}')

	@property
	def poses(self):
		"""Get the poses associated to this set of compounds"""
		from .pset import PoseSet
		ids = self.db.select_where(query='pose_id', table='pose', key=f'pose_compound in {tuple(self.ids)}', multiple=True)
		ids = [v for v, in ids]
		return PoseSet(self.db, ids)

	### FILTERING

	def get_by_tag(self,tag):
		"""Get all child compounds with a certain tag"""
		values = self.db.select_where(query='tag_compound', table='tag', key='name', value=tag, multiple=True)
		ids = [v for v, in values if v and v in self.ids]
		return CompoundSet(self.db, ids)

	def get_by_metadata(self, key: str, value: str | None = None):
		"""Get all child compounds with by their metadata. If no value is passed, then simply containing the key in the metadata dictionary is sufficient"""
		results = self.db.select(query='compound_id, compound_metadata', table='compound', multiple=True)
		if value is None:
			ids = [i for i,d in results if d and f'"{key}":' in d and i in self.ids]
		else:
			if isinstance(value, str):
				value = f'"{value}"'
			ids = [i for i,d in results if d and f'"{key}": {value}' in d and i in self.ids]
		return CompoundSet(self.db, ids)		

	### OUTPUT

	def draw(self):
		"""Draw a grid of all contained molecules"""
		from molparse.rdkit import draw_grid

		data = [(str(c), c.mol) for c in self]

		mols = [d[1] for d in data]
		labels = [d[0] for d in data]

		return draw_grid(mols, labels=labels)

	def summary(self):
		"""Print a summary of this compound set"""
		logger.header('CompoundSet()')
		logger.var('#compounds', len(self))
		logger.var('#poses', self.num_poses)
		logger.var('tags', self.tags)

	def interactive(self):
		"""Creates a ipywidget to interactively navigate this PoseSet."""

		from ipywidgets import interactive, BoundedIntText, Checkbox, interactive_output, HBox, GridBox, Layout, VBox
		from IPython.display import display
		from pprint import pprint

		a = BoundedIntText(
				value=0,
				min=0,
				max=len(self)-1,
				step=1,
				description=f'Comp (/{len(self)}):',
				disabled=False,
			)

		b = Checkbox(description='Name', value=True)
		c = Checkbox(description='Summary', value=False)
		d = Checkbox(description='2D', value=True)
		e = Checkbox(description='poses', value=False)
		f = Checkbox(description='Metadata', value=False)

		ui = GridBox([b, c, d, e, f], layout=Layout(grid_template_columns="repeat(5, 100px)"))
		ui = VBox([a, ui])
		
		def widget(i, name=True, summary=True, draw=True, poses=True, metadata=True):
			comp = self[i]
			
			if name and not summary:
				print(repr(comp))

			if summary: 
				comp.summary(metadata=False, draw=False)
			
			if draw: 
				comp.draw()
			
			if poses and comp.poses: 
				comp.poses.draw()
			
			if metadata:
				logger.title('Metadata:')
				pprint(comp.metadata)

		out = interactive_output(widget, {'i': a, 'name': b, 'summary': c, 'draw':d, 'poses':e, 'metadata':f})

		display(ui, out)


	### OTHER METHODS

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
		return self.db.get_compound(table=self.table, id=index)

class IngredientSet:

	"""Refactor as a wrapper to a dataframe???"""

	_columns = [
		'compound_id',
		'amount',
		'quote_id',
		'supplier',
		'max_lead_time',
	]
	
	def __init__(self, db, ingredients=None):

		from pandas import DataFrame

		ingredients = ingredients or []

		self._db = db
		
		self._data = DataFrame(columns=self._columns)
		
		for ingredient in ingredients:
			self.add(ingredient)

	### PROPERTIES

	@property
	def df(self):
		return self._data

	### METHODS

	def add(self, ingredient=None, *, compound_id=None, amount=None, quote_id=None, supplier=None, max_lead_time=None):

		from pandas import DataFrame, concat

		if ingredient:
			assert ingredient._table == 'ingredient'
			compound_id = ingredient.compound_id
			amount = ingredient.amount
			quote_id = ingredient.quote_id
			supplier = ingredient.supplier
			max_lead_time = ingredient.max_lead_time
			
		else:
			assert id
			assert amount

		self._data = concat([self._data, DataFrame([dict(compound_id=compound_id, amount=amount, quote_id=quote_id, supplier=supplier, max_lead_time=max_lead_time)])], ignore_index=True)

	def get_ingredient(self, series):
		return Ingredient(
			db=self._db,
			compound=series['compound_id'],
			amount=series['amount'],
			quote=series['quote_id'],
			supplier=series['supplier'],
			max_lead_time=series['max_lead_time'],
		)

	### DUNDERS

	def __len__(self):
		return len(self._data)

	def __repr__(self):
		return f'{mcol.bold}{mcol.underline}''{'f'I x {len(self)}''}'f'{mcol.unbold}{mcol.ununderline}'

	def __add__(self, other):
		for id, amount in other._data.items():
			self.add(id=id, amount=amount)
		return self

	def __getitem__(self, key):
		match key:
			case int():
				series = self.df.loc[key]
				self.get_ingredient(series)
		
			case _:
				raise NotImplementedError

	def __iter__(self):
		return iter(self.get_ingredient(s) for i,s in self.df.iterrows())