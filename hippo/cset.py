
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

	@property
	def reactants(self):
		"""Returns a CompoundSet of all compounds that are used as a reactants"""
		# ids = self.db.select(table='reactant', query='DISTINCT reactant_compound', multiple=True)
		ids = self.db.execute('SELECT reactant_compound FROM reactant LEFT JOIN reaction ON reactant.reactant_compound = reaction.reaction_product WHERE reaction.reaction_product IS NULL').fetchall()
		ids = [q for q, in ids]
		from .cset import CompoundSet
		return CompoundSet(self.db, ids)

	@property
	def products(self):
		"""Returns a CompoundSet of all compounds that are a product of a reaction but not a reactant"""
		ids = self.db.execute('SELECT reaction_product FROM reaction LEFT JOIN reactant ON reaction.reaction_product = reactant.reactant_compound WHERE reactant.reactant_compound IS NULL').fetchall()
		ids = [q for q, in ids]
		from .cset import CompoundSet
		return CompoundSet(self.db, ids)

	@property
	def intermediates(self):
		"""Returns a CompoundSet of all compounds that are products and reactants"""
		ids = self.db.execute('SELECT DISTINCT reaction_product FROM reaction INNER JOIN reactant ON reaction.reaction_product = reactant.reactant_compound').fetchall()
		ids = [q for q, in ids]
		from .cset import CompoundSet
		return CompoundSet(self.db, ids)

	@property
	def num_reactants(self):
		"""Returns the number of reactants (see HIPPO.reactants)"""
		return len(self.reactants)

	@property
	def num_intermediates(self):
		"""Returns the number of intermediates (see HIPPO.intermediates)"""
		return len(self.intermediates)
	
	@property
	def num_products(self):
		"""Returns the number of products (see HIPPO.products)"""
		return len(self.products)

	@property
	def elabs(self):
		"""Returns a CompoundSet of all compounds that are a an elaboration of an existing base"""
		ids = self.db.select_where(query='compound_id', table='compound', key='compound_base IS NOT NULL', multiple=True)
		ids = [q for q, in ids]
		from .cset import CompoundSet
		return CompoundSet(self.db, ids)

	@property
	def bases(self):
		"""Returns a CompoundSet of all compounds that are the basis for a set of elaborations"""
		ids = self.db.select_where(query='DISTINCT compound_base', table='compound', key='compound_base IS NOT NULL', multiple=True, none='quiet')
		ids = [q for q, in ids]
		from .cset import CompoundSet
		return CompoundSet(self.db, ids)

	@property
	def num_elabs(self):
		"""Returns the number of compounds that are a an elaboration of an existing base"""
		return len(self.elabs)

	@property
	def num_bases(self):
		"""Returns the number of compounds that are the basis for a set of elaborations"""
		return len(self.bases)

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
		logger.var('#bases', self.num_bases)
		logger.var('#elabs', self.num_elabs)
		logger.var('#reactants', self.num_reactants)
		logger.var('#intermediates', self.num_intermediates)
		logger.var('#products', self.num_products)

	def draw(self):
		return self[self.ids].draw()

	def interactive(self):
		return self[self.ids].interactive()

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
		return f'{mcol.bold}{mcol.underline}''{'f'C x {len(self)}''}'f'{mcol.unbold}{mcol.ununderline}'

	def __len__(self) -> int:
		return self.db.count(self.table)

	def __iter__(self):
		return iter(self[i+1] for i in range(len(self)))


class CompoundSet:
	"""Object representing a subset of the 'compound' table in the :class:`.Database`."""

	_table = 'compound'
	
	def __init__(self,
		db: Database,
		indices: list = None,
	):

		self._db = db

		indices = indices or []

		if not isinstance(indices, list):
			indices = list(indices)

		assert all(isinstance(i, int) for i in indices)

		self._indices = sorted(list(set(indices)))

	### PROPERTIES

	@property
	def db(self):
		return self._db

	@property
	def table(self):
		return self._table

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
		result = self.db.select_where(query='compound_alias', table='compound', key=f'compound_id in {self.str_ids}', multiple=True)
		return [q for q, in result]

	@property
	def smiles(self) -> list:
		"""Returns the smiles of child compounds"""
		result = self.db.select_where(query='compound_smiles', table='compound', key=f'compound_id in {self.str_ids}', multiple=True)
		return [q for q, in result]
	
	@property
	def inchikeys(self):
		"""Returns the inchikeys of compounds in this set"""
		result = self.db.select_where(query='compound_inchikey', table='compound', key=f'compound_id in {self.str_ids}', multiple=True)
		return [q for q, in result]

	@property
	def tags(self):
		"""Returns the set of unique tags present in this compound set"""
		values = self.db.select_where(table='tag', query='DISTINCT tag_name', key=f'tag_compound in {self.str_ids}', multiple=True)
		if not values:
			return set()
		return set(v for v, in values)

	@property
	def num_poses(self):
		"""Count the poses associated to this set of compounds"""
		from .pset import PoseSet
		return self.db.count_where(table='pose', key=f'pose_compound in {self.str_ids}')

	@property
	def poses(self):
		"""Get the poses associated to this set of compounds"""
		from .pset import PoseSet
		ids = self.db.select_where(query='pose_id', table='pose', key=f'pose_compound in {self.str_ids}', multiple=True)
		ids = [v for v, in ids]
		return PoseSet(self.db, ids)

	@property
	def best_placed_poses(self):
		"""Get the best placed pose for each compound in this set"""
		from .pset import PoseSet
		query = self.db.select_where(table='pose', query='pose_id, MIN(pose_distance_score)', key=f'pose_compound in {self.str_ids} GROUP BY pose_compound', multiple=True)
		ids = [i for i,s in query]
		return PoseSet(self.db, ids)

	@property
	def str_ids(self):
		return str(tuple(self.ids)).replace(',)',')')

	@property
	def num_heavy_atoms(self):
		return sum([c.num_heavy_atoms for c in self])

	@property
	def num_rings(self):
		return sum([c.num_rings for c in self])

	@property
	def formula(self):
		from .tools import atomtype_dict_to_formula
		return atomtype_dict_to_formula(self.atomtype_dict)

	@property
	def atomtype_dict(self):
		from .tools import formula_to_atomtype_dict, combine_atomtype_dicts
		atomtype_dicts = [c.atomtype_dict for c in self]
		return combine_atomtype_dicts(atomtype_dicts)

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

	def add(self, compound):
		"""Add a compound to this set"""
		if isinstance(compound, Compound):
			compound = compound.id

		if compound not in self.ids:
			from bisect import insort
			insort(self.ids, compound)
	
	def get_recipe(self, 
		amount: float = 1, 
		simplify_network: bool = True,
	):
		"""Generate the :class:`.Recipe` to make these compounds."""

		n_comps = len(self)

		assert simplify_network

		if not hasattr(amount, '__iter__'):
			amount = [amount] * n_comps

		from .cset import IngredientSet
		from .rset import ReactionSet
		
		products = IngredientSet(self.db)
		reactants = IngredientSet(self.db)
		intermediates = IngredientSet(self.db)
		reactions = ReactionSet(self.db)

		for comp, a in zip(self, amount):

			reax = comp.reactions

			assert len(reax) == 1
			# logger.warning(f'Picking first reaction for {comp} from {reax}')

			recipe = reax[0].get_recipe(a)

			# recipe.summary()
			
			products.add(comp.as_ingredient(amount=a))

			intermediates += recipe.intermediates
			
			reactants += recipe.reactants

			reactions += recipe.reactions
			
		return Recipe(products=products, reactants=reactants, reactions=reactions, intermediates=intermediates)

	def copy(self):
		return CompoundSet(self.db, self.ids)
	
	def get_df(self, **kwargs):

		from tqdm import tqdm
		from pandas import DataFrame

		data = []

		for comp in tqdm(self):
			d = comp.get_dict(**kwargs)
			data.append(d)

		return DataFrame(data)

	### DUNDERS

	# def __repr__(self) -> str:
		# return f'{mcol.bold}{mcol.underline}''{'f'C x {len(self)}''}'f'{mcol.unbold}{mcol.ununderline}'

	def __len__(self) -> int:
		return len(self.indices)

	def __iter__(self):
		return iter(self.db.get_compound(table=self.table, id=i) for i in self.indices)

	def __getitem__(self, key) -> Compound:
		match key:
			case int():
				index = self.indices[key]
				return self.db.get_compound(table=self.table, id=index)
			
			case slice():
				indices = self.indices[key]
				return	CompoundSet(self.db, indices)

			case _:
				raise NotImplementedError	

	def __sub__(self, other):

		match other:
		
			case CompoundSet():
				ids = set(self.ids) - set(other.ids)
				return CompoundSet(self.db, ids)
				
			case _:
				raise NotImplementedError

	def __add__(self, other):

		match other:

			case Compound():
				return self.add(other)

			case int():
				return self.add(other)
		
			case CompoundSet():
				ids = set(self.ids) | set(other.ids)
				return CompoundSet(self.db, ids)
				
			case _:
				raise NotImplementedError

	def __xor__(self, other):

		match other:
		
			case CompoundSet():
				ids = set(self.ids) ^ set(other.ids)
				return CompoundSet(self.db, ids)
				
			case _:
				raise NotImplementedError

	def __repr__(self) -> str:
		return f'{mcol.bold}{mcol.underline}''{'f'C x {len(self)}''}'f'{mcol.unbold}{mcol.ununderline}'


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

	@property
	def db(self):
		return self._db

	### METHODS

	def add(self, ingredient=None, *, compound_id=None, amount=None, quote_id=None, supplier=None, max_lead_time=None):

		from pandas import DataFrame, concat

		if ingredient:
			assert ingredient._table == 'ingredient'
			compound_id = ingredient.compound_id
			amount = ingredient.amount

			if ingredient.quote and not ingredient.quote_id:
				logger.warning(f'Losing quote! {ingredient.quote=}')

			quote_id = ingredient.quote_id
			
			supplier = ingredient.supplier
			max_lead_time = ingredient.max_lead_time
			
		else:
			assert compound_id
			assert amount

		if self._data.empty:
			addition = DataFrame([dict(compound_id=compound_id, amount=amount, quote_id=quote_id, supplier=supplier, max_lead_time=max_lead_time)])
			self._data = addition

		else:

			if compound_id in self._data['compound_id'].values:
				index = self._data.index[self._data['compound_id'] == compound_id].tolist()[0]
				self._data.loc[index, 'amount'] += amount
				self._data.loc[index, 'quote_id'] = None

			else:
				from numpy import nan
				addition = DataFrame([dict(compound_id=compound_id, amount=amount, quote_id=quote_id, supplier=supplier, max_lead_time=max_lead_time)])
				self._data = concat([self._data, addition], ignore_index=True)
				self._data = self._data.replace({nan: None})

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
		for i,row in other._data.iterrows():
			self.add(
				compound_id=row.compound_id,
				amount=row.amount,
				quote_id=row.quote_id,
				supplier=row.supplier,
				max_lead_time=row.max_lead_time,
			)
		return self

	def __getitem__(self, key):
		match key:
			case int():
				series = self.df.loc[key]
				return self.get_ingredient(series)
		
			case _:
				raise NotImplementedError

	def __iter__(self):
		return iter(self.get_ingredient(s) for i,s in self.df.iterrows())

	def __call__(self, *, compound_id):

		if compound_id:

			# get the ingredient with the matching compound ID
			matches = self.df[self.df['compound_id'] == compound_id]

			if len(matches) != 1:

				logger.warning(f'Multiple ingredients in set with {compound_id=}')
				# print(matches)

				return IngredientSet(self.db, [self.get_ingredient(s) for i,s in matches.iterrows()])

			return self.get_ingredient(matches.iloc[0])

		else:
			raise NotImplementedError
