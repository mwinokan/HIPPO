
# from .tools import df_row_to_dict

from .compound import Compound, Ingredient
from .db import Database

from .recipe import Recipe

from numpy import int64, nan, isnan

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
		""" """
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
		"""Get all child compounds with a certain tag

		:param tag: 

		"""
		values = self.db.select_where(query='tag_compound', table='tag', key='name', value=tag, multiple=True)
		ids = [v for v, in values if v]
		return self[ids]

	def get_by_metadata(self, key: str, value: str | None = None):
		"""Get all child compounds with by their metadata. If no value is passed, then simply containing the key in the metadata dictionary is sufficient

		:param key: str:
		:param value: str | None:  (Default value = None)
		:param key: str: 
		:param value: str | None:  (Default value = None)

		"""
		results = self.db.select(query='compound_id, compound_metadata', table='compound', multiple=True)
		if value is None:
			ids = [i for i,d in results if d and f'"{key}":' in d]
		else:
			if isinstance(value, str):
				value = f'"{value}"'
			ids = [i for i,d in results if d and f'"{key}": {value}' in d]
		return self[ids]		

	def get_by_base(self,base):
		"""

		:param base: 

		"""

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
		""" """
		return self[self.ids].draw()

	def interactive(self):
		""" """
		return self[self.ids].interactive()

	### DUNDERS

	def __call__(self, *, tag=None, base=None):
		if tag:
			return self.get_by_tag(tag)
		elif base:
			return self.get_by_base(base)
		else:
			raise NotImplementedError(f'{type(i)=}')

	def __getitem__(self, key) -> Compound:
		
		match key:

			# case int():
			case key if isinstance(key, int) or isinstance(key, int64):

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

			case key if isinstance(key, list) or isinstance(key, tuple) or isinstance(key, set):

				indices = []
				for i in key:
					if isinstance(i,int) or isinstance(i, int64):
						index = i
					elif isinstance(i,str):
						index = self.db.get_compound_id(name=i)
					else:
						raise NotImplementedError

					assert index
					indices.append(index)

				return CompoundSet(self.db, indices)

			case slice():
				ids = self.db.slice_ids(table=self.table, start=key.start, stop=key.stop, step=key.step)
				return self[ids]

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
		sort: bool = True,
	):

		self._db = db

		indices = indices or []

		if not isinstance(indices, list):
			indices = list(indices)

		indices = [int(i) for i in indices]

		if sort:
			self._indices = sorted(list(set(indices)))
		else:
			self._indices = list(set(indices))

	### PROPERTIES

	@property
	def db(self):
		""" """
		return self._db

	@property
	def table(self):
		""" """
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
		""" """
		return str(tuple(self.ids)).replace(',)',')')

	@property
	def num_heavy_atoms(self):
		""" """
		return sum([c.num_heavy_atoms for c in self])

	@property
	def num_rings(self):
		""" """
		return sum([c.num_rings for c in self])

	@property
	def formula(self):
		""" """
		from .tools import atomtype_dict_to_formula
		return atomtype_dict_to_formula(self.atomtype_dict)

	@property
	def atomtype_dict(self):
		""" """
		from .tools import formula_to_atomtype_dict, combine_atomtype_dicts
		atomtype_dicts = [c.atomtype_dict for c in self]
		return combine_atomtype_dicts(atomtype_dicts)

	@property
	def num_atoms_added(self) -> list[int]:

		"""Calculate the number of atoms added w.r.t the base

		:returns: list of number of atoms added values

		"""

		query = self.db.execute(f"""
		WITH nums AS (
			SELECT A.compound_id AS comp_id, 
			mol_num_hvyatms(A.compound_mol) - mol_num_hvyatms(B.compound_mol) AS diff 
			FROM compound A, compound B
			WHERE A.compound_base = B.compound_id
			AND A.compound_id IN {self.str_ids}
		)

		SELECT compound_id, diff FROM compound
		LEFT JOIN nums
		ON comp_id = compound_id
		WHERE compound_id IN {self.str_ids}
		""").fetchall()

		lookup = {k:v for k,v in query}

		return [lookup[i] for i in self.indices]

	### FILTERING

	def get_by_tag(self,tag):
		"""Get all child compounds with a certain tag

		:param tag: 

		"""
		values = self.db.select_where(query='tag_compound', table='tag', key='name', value=tag, multiple=True)
		ids = [v for v, in values if v and v in self.ids]
		return CompoundSet(self.db, ids)

	def get_by_metadata(self, key: str, value: str | None = None):
		"""Get all child compounds with by their metadata. If no value is passed, then simply containing the key in the metadata dictionary is sufficient

		:param key: str:
		:param value: str | None:  (Default value = None)
		:param key: str: 
		:param value: str | None:  (Default value = None)

		"""
		results = self.db.select(query='compound_id, compound_metadata', table='compound', multiple=True)
		if value is None:
			ids = [i for i,d in results if d and f'"{key}":' in d and i in self.ids]
		else:
			if isinstance(value, str):
				value = f'"{value}"'
			ids = [i for i,d in results if d and f'"{key}": {value}' in d and i in self.ids]
		return CompoundSet(self.db, ids)		

	def get_all_possible_reactants(self, debug=False):
		"""Recursively searches for all the reactants that could possible be needed to synthesise these compounds.

		:param debug: Default value = False)

		"""
		all_reactants, all_reactions = self.db.get_unsolved_reaction_tree(product_ids=self.ids, debug=debug) 
		return all_reactants

	def get_all_possible_reactions(self, debug=False):
		"""Recursively searches for all the reactants that could possible be needed to synthesise these compounds.

		:param debug: Default value = False)

		"""
		all_reactants, all_reactions = self.db.get_unsolved_reaction_tree(product_ids=self.ids, debug=debug) 
		return all_reactions

	### CONSOLE / NOTEBOOK OUTPUT

	def draw(self):
		"""Draw a grid of all contained molecules"""
		from molparse.rdkit import draw_grid

		data = [(str(c), c.mol) for c in self]

		mols = [d[1] for d in data]
		labels = [d[0] for d in data]

		display(draw_grid(mols, labels=labels))

	def grid(self):
		""" """
		self.draw()

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
		f = Checkbox(description='reactions', value=False)
		g = Checkbox(description='Metadata', value=False)

		ui = GridBox([b, c, d, e, f, g], layout=Layout(grid_template_columns="repeat(6, 100px)"))
		ui = VBox([a, ui])
		
		def widget(i, name=True, summary=True, draw=True, poses=True, reactions=True, metadata=True):
			"""

			:param i: param name:  (Default value = True)
			:param summary: Default value = True)
			:param draw: Default value = True)
			:param poses: Default value = True)
			:param reactions: Default value = True)
			:param metadata: Default value = True)
			:param name:  (Default value = True)

			"""
			comp = self[i]
			
			if name and not summary:
				print(repr(comp))

			if summary: 
				comp.summary(metadata=False, draw=False)
			
			if draw: 
				comp.draw()
			
			if poses and (pset := comp.poses): 
				for p in pset:
					print(repr(p))
				pset.draw()

			if reactions and (reactions := comp.reactions):
				for r in reactions:
					print(repr(r))
					r.draw()
			
			if metadata:
				logger.title('Metadata:')
				pprint(comp.metadata)

		out = interactive_output(widget, {'i': a, 'name': b, 'summary': c, 'draw':d, 'poses':e, 'reactions':f, 'metadata':g})

		display(ui, out)		

	### OTHER METHODS

	def add(self, compound):
		"""Add a compound to this set

		:param compound: 

		"""
		if isinstance(compound, Compound):
			compound = compound.id

		if compound not in self.ids:
			from bisect import insort
			insort(self.ids, compound)
	
	def get_recipes(self, 
		amount: float = 1, 
		debug=False,
		pick_cheapest: bool = False,
		permitted_reactions=None,
		quoted_only: bool = False,
		supplier: None | str = None,
		**kwargs,
	):
		"""Generate the :class:`.Recipe` to make these compounds.

		:param amount: float:  (Default value = 1)
		:param debug: Default value = False)
		:param pick_cheapest: bool:  (Default value = False)
		:param permitted_reactions: Default value = None)
		:param quoted_only: bool:  (Default value = False)
		:param supplier: None | str:  (Default value = None)
		:param amount: float:  (Default value = 1)
		:param pick_cheapest: bool:  (Default value = False)
		:param quoted_only: bool:  (Default value = False)
		:param supplier: None | str:  (Default value = None)

		"""
		from .recipe import Recipe
		return Recipe.from_compounds(self, 
			amount=amount, 
			debug=debug, 
			pick_cheapest=pick_cheapest, 
			permitted_reactions=permitted_reactions, 
			quoted_only=quoted_only, 
			supplier=supplier, 
			**kwargs
		)

	def copy(self):
		"""Returns a copy of this set"""
		return CompoundSet(self.db, self.ids)

	def shuffled(self):
		"""Returns a randomised copy of this set"""
		copy = self.copy()
		copy.shuffle()
		return copy

	def pop(self):
		"""Pop the last compound in this set"""
		c_id = self.pop_id()
		return self.db.get_compound(id=c_id)

	def pop_id(self):
		"""Pop the last compound id in this set"""
		return self._indices.pop()

	def shuffle(self):
		"""Randomises the order of compounds in this set"""
		from random import shuffle
		shuffle(self._indices)
	
	def get_df(self, mol=False, reactions=False, metadata=False, count_by_target=False, poses=False, **kwargs):
		"""

		:param mol: Default value = False)
		:param reactions: Default value = False)
		:param metadata: Default value = False)
		:param count_by_target: Default value = False)
		:param poses: Default value = False)

		"""

		from tqdm import tqdm
		from pandas import DataFrame

		data = []

		for comp in tqdm(self):
			d = comp.get_dict(mol=mol, reactions=reactions, metadata=metadata, count_by_target=count_by_target, poses=poses, **kwargs)
			data.append(d)

		return DataFrame(data)

	def get_quoted(self, *, supplier='any'):
		"""Get all member compounds that have a quote from given supplier

		:param *: 
		:param supplier:  (Default value = 'any')

		"""

		if supplier != 'any':
			key = f'quote_compound IN {self.str_ids} AND quote_supplier = "{supplier}"'
		else:
			key = f'quote_compound IN {self.str_ids}'

		ids = self.db.select_where(
			table='quote', 
			query='DISTINCT quote_compound', 
			key=key, 
			multiple=True,
		)
		
		ids = [i for i, in ids]
		return CompoundSet(self.db, ids)

	def get_unquoted(self, *, supplier='any'):
		"""Get all member compounds that do not have a quote from given supplier

		:param *: 
		:param supplier:  (Default value = 'any')

		"""
		quoted = self.get_quoted(supplier=supplier)
		return self - quoted

	def get_dict(self):
		"""Get a dictionary object with all serialisable data needed to reconstruct this set"""
		return dict(db=str(self.db.path.resolve()), indices=self.indices)

	def write_smiles_csv(self, file):
		"""

		:param file: 

		"""
		from pandas import DataFrame
		smiles = self.smiles
		df = DataFrame(dict(smiles=smiles))
		df.to_csv(file)

	def write_postera_csv(self, file, *, supplier='Enamine', prefix='fragment'):
		"""

		:param file:
		:param supplier: Default value = 'Enamine')
		:param prefix: Default value = 'fragment')

		"""
		from datetime import date as dt
		from pandas import DataFrame
		from tqdm import tqdm

		if prefix:
			prefix = f'{prefix}_'

		data = []

		for c in tqdm(self, total=len(self)):
			
			# get props
			smiles = c.smiles
			tags = c.tags
			metadata = c.metadata
			poses = c.poses
			base = c.base

			# method
			assert len(tags) == 1, c
			method = tags[0]

			# date
			date = dt.today()

			# author
			assert 'author' in metadata, c
			author = metadata['author']


			match len(poses):
				case 1:
					pose = poses[0]
				case 0:
					logger.warning(f'{c} has no poses')
					assert base
					pose = base.poses[0]
				case _:
					logger.warning(f'{c} has multiple poses')
					pose = poses[0]

			# extract inspirations
			inspirations = pose.inspirations
			inspiration_names = ','.join(inspirations.names)
			inspiration_smiles = '.'.join(inspirations.smiles)

			# quote info
			quotes = c.get_quotes(supplier=supplier)
			assert len(quotes) == 1, c
			quote = quotes[0]
			catalog_id = quote.entry
			catalog_price = quote.price
			catalog_lead_time = quote.lead_time

			# hippo string
			hippo_str = f'compound={c.id}, pose={pose.id}'

			# create row
			data.append({
				"SMILES":smiles,
				f"{prefix}HIPPO_IDs":hippo_str,
				f"{prefix}method":method,
				f"{prefix}export_date":date,
				f"{prefix}author":author,
				f"{prefix}inspiration_names":inspiration_names,
				f"{prefix}inspiration_SMILES":inspiration_smiles,
				f"{prefix}supplier":supplier,
				f"{prefix}supplier_catalogue":quote.catalogue,
				f"{prefix}supplier_ID":catalog_id,
				f"{prefix}supplier_price":catalog_price,
				f"{prefix}supplier_lead_time":catalog_lead_time,
			})

		df = DataFrame(data)

		logger.writing(file)
		df.to_csv(file, index=False)

		return df

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

			case IngredientSet():
				logger.warning('Subtracting IngredientSet from CompoundSet. Ignoring quote/amount data')
				ids = set(self.ids) - set([int(i) for i in other.compound_ids])
				# print(self.ids)
				# print(ids)
				# print([int(i) for i in other.compound_ids])
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

			case IngredientSet():
				ids = set(self.ids) | set(other.compound_ids)
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

	def __contains__(self, other):
		match other:
			case Compound():
				id = other.id
			case Ingredient():
				id = other.compound_id
			case int():
				id = other

		return id in set(self.ids)


class IngredientSet:

	"""Refactor as a wrapper to a dataframe???"""

	_columns = [
		'compound_id',
		'amount',
		'quote_id',
		'supplier',
		'max_lead_time',
		'quoted_amount',
	]
	
	def __init__(self, db, ingredients=None, supplier=None):

		from pandas import DataFrame

		ingredients = ingredients or []

		self._db = db
		
		self._data = DataFrame(columns=self._columns, dtype=object)

		self._supplier = supplier
		
		for ingredient in ingredients:
			self.add(ingredient)

	@classmethod
	def from_ingredient_df(cls, db, df, supplier=None):
		"""

		:param db: param df:
		:param supplier: Default value = None)
		:param df: 

		"""
		# from numpy import nan
		self = cls.__new__(cls)

		for col in cls._columns:
			if col not in df.columns:
				logger.debug(f'Adding column {col}')
				df[col] = None				
		
		self._db = db
		self._data = df.copy()
		self._supplier = supplier

		return self

	@classmethod
	def from_json(cls, db, path, supplier, data=None):
		"""

		:param db: param path:
		:param supplier: param data:  (Default value = None)
		:param path: 
		:param data:  (Default value = None)

		"""

		if not data:
			import json
			data = json.load(open(path,'rt'))

		from pandas import DataFrame
		df = DataFrame(columns=cls._columns, dtype=object)

		for col in cls._columns:
			df[col] = data[col]
		
		return cls.from_ingredient_df(db=db, df=df, supplier=supplier)

	@classmethod
	def from_ingredient_dicts(cls, db, dicts, supplier=None):
		"""

		:param db: param dicts:
		:param supplier: Default value = None)
		:param dicts: 

		"""
		from pandas import DataFrame
		df = DataFrame(dicts, dtype=object)
		# for col in cls._columns:
		# 	if col not in df.columns:
		# 		df[col] = None
		return cls.from_ingredient_df(db=db, df=df, supplier=supplier)

	@classmethod
	def from_compounds(cls, compounds, ids=None, db=None, amount=1):
		"""

		:param compounds: param ids:  (Default value = None)
		:param db: Default value = None)
		:param amount: Default value = 1)
		:param ids:  (Default value = None)

		"""

		from pandas import DataFrame

		if not ids:
			ids = compounds.ids

		if not db:
			db = compounds.db

		df = DataFrame(dict(compound_id=ids, amount=amount, quote_id=None, supplier=None, max_lead_time=None, quoted_amount=None), dtype=object)

		return cls.from_ingredient_df(db, df)

	### PROPERTIES

	@property
	def df(self):
		""" """
		return self._data

	@property
	def db(self):
		""" """
		return self._db

	@property
	def price_df(self):
		""" """
		df = self.df.copy()
		df['price'] = [i.price for i in self]
		return df

	@property
	def price(self):
		""" """
		return self.get_price()

	@property
	def supplier(self):
		""" """
		return self._supplier

	@supplier.setter
	def supplier(self, s):
		"""

		:param s: 

		"""
		
		if isinstance(s, list) or isinstance(s, tuple):
			for x in s:
				assert isinstance(x, str)
		else:
			assert isinstance(s, str)

		self._supplier = s
		self.df['supplier'] = [s] * len(self)

	@property
	def smiles(self):
		""" """
		compound_ids = list(self.df['compound_id'])
		result = self.db.select_where(query='compound_smiles', table='compound', key=f'compound_id in {tuple(compound_ids)}', multiple=True)
		return [q for q, in result]

	@property
	def inchikeys(self):
		""" """
		compound_ids = list(self.df['compound_id'])
		result = self.db.select_where(query='compound_inchikey', table='compound', key=f'compound_id in {tuple(compound_ids)}', multiple=True)
		return [q for q, in result]

	@property
	def compound_ids(self):
		""" """
		return list(self.df['compound_id'].values)

	@property
	def ids(self):
		""" """
		return self.compound_ids

	@property
	def str_compound_ids(self):
		""" """
		return str(tuple(self.df['compound_id'].values)).replace(',)',')')

	@property
	def compounds(self):
		""" """
		return CompoundSet(self.db, self.compound_ids)

	### METHODS

	def get_price(self, supplier=None):
		"""

		:param supplier: Default value = None)

		"""

		from .price import Price

		pairs = { i:q for i,q in enumerate(self.df['quote_id'])}

		quote_ids = [q for q in pairs.values() if q is not None and not isnan(q)]

		if quote_ids:

			quote_id_str = str(tuple(quote_ids)).replace(',)', ')')

			if supplier:
				result = self.db.select_where(query='quote_price, quote_currency', table='quote', key=f'quote_id in {quote_id_str} AND quote_supplier = "{supplier}"', multiple=True)
			else:
				result = self.db.select_where(query='quote_price, quote_currency', table='quote', key=f'quote_id in {quote_id_str}', multiple=True)

			prices = [Price(a,b) for a,b in result]
			quoted = sum(prices, Price.null())
			
		else:

			quoted = Price.null()

		unquoted = [i for i,q in pairs.items() if q is None or isnan(q)]

		unquoted_price = Price.null()
		for i in unquoted:
			ingredient = self[i]
			p = ingredient.price
			unquoted_price += p
			
			quote = ingredient.quote
			self.df.loc[i, 'quote_id'] = quote.id
			assert quote.amount
			self.df.loc[i, 'quoted_amount'] = quote.amount
		
		return quoted + unquoted_price

	def interactive(self, **kwargs):
		"""Wrapper for :meth:`.CompoundSet.interactive`"""
		return self.compounds.interactive(**kwargs)

	def add(self, ingredient=None, *, compound_id=None, amount=None, quote_id=None, supplier=None, max_lead_time=None, quoted_amount=None, debug=False):
		"""

		:param ingredient: Default value = None)
		:param *: 
		:param compound_id:  (Default value = None)
		:param amount:  (Default value = None)
		:param quote_id:  (Default value = None)
		:param supplier:  (Default value = None)
		:param max_lead_time:  (Default value = None)
		:param quoted_amount:  (Default value = None)
		:param debug:  (Default value = False)

		"""

		from pandas import DataFrame, concat

		if ingredient:
			assert ingredient._table == 'ingredient'
			compound_id = ingredient.compound_id
			amount = ingredient.amount

			if (q := ingredient.quote) and not ingredient.quote_id:
				logger.warning(f'Losing quote! {ingredient.quote=}')
			
			supplier = ingredient.supplier
			max_lead_time = ingredient.max_lead_time

			quote_id = q.id
			quoted_amount = q.amount
			
		else:
			assert compound_id
			assert amount

		if quote_id:
			assert quoted_amount

		# if self.supplier:
			# assert supplier == self.supplier, f'current: {self.supplier=}, adding: {supplier=}'

		supplier = self.supplier

		if self._data.empty:
			addition = DataFrame([dict(compound_id=compound_id, amount=amount, quote_id=quote_id, supplier=supplier, max_lead_time=max_lead_time)], dtype=object)
			self._data = addition
			# self._data = self._data.replace({nan: None})

		else:

			if compound_id in self._data['compound_id'].values:
				index = self._data.index[self._data['compound_id'] == compound_id].tolist()[0]
				self._data.loc[index, 'amount'] += amount

				# discard if the quote is no longer valid
				if (a := self.df.loc[index, 'quoted_amount']) and a < self.df.loc[index, 'amount']:
					self._data.loc[index, 'quote_id'] = None
					self._data.loc[index, 'quoted_amount'] = None

				if debug and supplier:
					logger.debug('Adding to existing ingredient')
					logger.debug(f'{self._data.loc[index, "supplier"]=}')
					logger.debug(f'{supplier=}')

			else:
				# from numpy import nan
				addition = DataFrame([dict(compound_id=compound_id, amount=amount, quote_id=quote_id, supplier=supplier, max_lead_time=max_lead_time, quoted_amount=quoted_amount)], dtype=object)
				# print(self._data)
				# print(addition)
				# print(self._data.dtypes)
				# print(addition.dtypes)

				self._data = concat([self._data, addition], ignore_index=True, join='inner')
				# self._data = self._data.replace({nan: None})

				if debug:
					logger.out(addition)

	def get_ingredient(self, series):
		"""

		:param series: 

		"""
		
		q_id = series['quote_id']

		if isinstance(q_id, float) and isnan(q_id):
			q_id = None

		return Ingredient(
			db=self._db,
			compound=series['compound_id'],
			amount=series['amount'],
			quote=q_id,
			supplier=series['supplier'],
			max_lead_time=series['max_lead_time'],
		)

	def copy(self):
		""" """
		return IngredientSet.from_ingredient_df(self.db, self.df, supplier=self.supplier)

	def draw(self):
		""" """
		self.compounds.draw()

	def set_amounts(self, amount):
		"""

		:param amount: 

		"""

		self.df['amount'] = amount

		# if amounts are modified the quotes should be cleared
		self.df['quote_id'] = None

		assert all(self.df['supplier'].isna()) and all(self.df['max_lead_time'].isna())

		# update quotes
		pairs = self.db.execute(f'''
			WITH matching_quotes AS (
				SELECT quote_id, quote_compound, MIN(quote_price) FROM quote
				WHERE quote_compound IN {self.str_compound_ids}
				AND quote_amount >= {amount}
				GROUP BY quote_compound
			)
			SELECT compound_id, quote_id FROM compound
			LEFT JOIN matching_quotes ON quote_compound = compound_id
			WHERE compound_id IN {self.str_compound_ids}
		''').fetchall()

		for compound_id, quote_id in pairs:
			match = self.df.index[self.df['compound_id'] == compound_id][0]
			self.df.loc[match, 'quote_id'] = quote_id

	def get_dict(self, data_orient='list'):
		"""Get serialisable dictionary

		:param data_orient: Default value = 'list')

		"""
		return dict(db=str(self.db), supplier=self.supplier, data=self.df.to_dict(orient=data_orient))

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
		# self._data = self._data.replace({nan: None})
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
