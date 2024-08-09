
from .recipe import Recipe
from .cset import CompoundSet, IngredientSet

from tqdm import tqdm

from pathlib import Path
import json

import logging
logger = logging.getLogger('HIPPO')

class RandomRecipeGenerator:

	def __init__(self, db, *, 
		max_lead_time = None,
		max_reactions = None,
		suppliers: list | None = None,
		start_with: Recipe | CompoundSet | IngredientSet = None,
	):

		logger.debug('RandomRecipeGenerator.__init__()')

		# Static parameters
		self._db_path = db.path
		self._max_lead_time = max_lead_time
		self._suppliers = suppliers
		self._starting_recipe = start_with

		logger.var('database', self.db_path)
		logger.var('max_lead_time', self.max_lead_time)
		logger.var('suppliers', self.suppliers)

		# Database set up
		self._db = db

		# JSON I/O set up
		self._data_path = Path(str(self.db_path).replace('.sqlite','_rgen.json'))
		if self.data_path.exists():
			logger.warning(f'Will overwrite existing rgen data file: {self.data_path}')

		# caches
		self._reaction_checking_cache = {}
		self._reaction_reactant_cache = {}

		# Reactant pool
		self._reactant_pool = self.get_reactant_pool()
		logger.var('reactant_pool', self.reactant_pool)

		# Product pool
		logger.debug('Solving product pool...')
		self._product_pool = self.get_product_pool()

		# Remove starting recipe products
		if start_with:
			self._product_pool -= start_with.products.compounds
		
		logger.var('product_pool', self.product_pool)

		# Route pool
		logger.debug('Solving route pool...')
		self._route_pool = self.get_route_pool()

		# dump data
		self.dump_data()

	### FACTORIES

	@classmethod
	def from_json(cls, db, path):

		data = json.load(open(path, 'rt'))

		self = cls.__new__(cls)

		self._db_path = data['db_path']
		self._max_lead_time = data['max_lead_time']
		self._suppliers = data['suppliers']

		self._starting_recipe = Recipe.from_json(db=db, path=None, data=data['starting_recipe'])

		logger.var('database', self.db_path)
		logger.var('max_lead_time', self.max_lead_time)
		logger.var('suppliers', self.suppliers)

		self._db = db

		# JSON I/O set up
		self._data_path = Path(path)

		# caches
		self._reaction_checking_cache = {}
		self._reaction_reactant_cache = {}

		# Reactant pool
		self._reactant_pool = CompoundSet(db=db, indices=data['reactant_pool']['indices'], sort=False)
		logger.var('reactant_pool', self.reactant_pool)
		
		# Product pool
		self._product_pool = CompoundSet(db=db, indices=data['product_pool']['indices'], sort=False)
		logger.var('product_pool', self.product_pool)
				
		# Route pool
		from .recipe import RouteSet
		self._route_pool = RouteSet.from_json(path=None, data=data['route_pool'], db=db)
		
		return self

	### PROPERTIES

	@property
	def starting_recipe(self):
		return self._starting_recipe

	@property
	def db(self):
		return self._db

	@property
	def db_path(self):
		return self._db_path

	@property
	def suppliers_str(self):
		return str(tuple(self.suppliers)).replace(',)',')')
	
	@property
	def suppliers(self):
		return self._suppliers
	
	@property
	def max_lead_time(self):
		return self._max_lead_time
	
	@property
	def reactant_pool(self):
		return self._reactant_pool

	@property
	def product_pool(self):
		return self._product_pool

	@property
	def route_pool(self):
		return self._route_pool
		
	@property
	def data_path(self):
		return self._data_path

	### POOL METHODS

	def get_reactant_pool(self):
		"""Get all purchaseable reactants with acceptable supplier and leadtime"""

		if self.max_lead_time:

			raise NotImplementedError

		elif self.suppliers:

			sql = f'''
			SELECT DISTINCT quote_compound FROM quote
			WHERE quote_supplier IN {self.suppliers_str} 
			ORDER BY quote_compound
			'''

		else:

			sql = f'''
			SELECT DISTINCT quote_compound FROM quote 
			ORDER BY quote_compound
			'''

		compound_ids = self.db.execute(sql).fetchall()
		compound_ids = [c for c, in compound_ids]

		return CompoundSet(self.db, compound_ids)

	def get_product_pool(self):
		"""Get all quoted products that can be made with the reactant pool"""
		reactants = self.get_reactant_pool()
		products = Recipe.from_reactants(reactants=reactants, debug=False, return_products=True)
		return products

	def get_route_pool(self, mini_test=False):

		route_ids = self.db.select_where(table='route', query='route_id', key=f'route_product IN {self.product_pool.str_ids}', multiple=True)

		if mini_test:
			route_ids = route_ids[:100]

		routes = [self.db.get_route(id=route_id) for route_id, in tqdm(route_ids)]

		from .recipe import RouteSet
		return RouteSet(self.db, routes)

	### FILE I/O METHODS

	def dump_data(self):

		data = {}

		data['db_path'] = str(self.db_path.resolve())
		data['max_lead_time'] = self.max_lead_time
		data['suppliers'] = self.suppliers
		# data['max_reactions'] = self.max_reactions
		data['starting_recipe'] = self.starting_recipe.get_dict(serialise_price=True)

		data['reactant_pool'] = self.reactant_pool.get_dict()
		data['product_pool'] = self.product_pool.get_dict()
		data['route_pool'] = self.route_pool.get_dict()

		logger.writing(self.data_path)
		json.dump(data, open(self.data_path, 'wt'), indent=4)

		# logger.debug('Clearing route_pool database pointers... ')
		# self.route_pool.clear_db_pointers()
		# logger.writing(self.data_path)
		# pickle.dump(self.route_pool, open(self.data_path, 'wb'))
		# self.route_pool.clear_db_pointers()
		# logger.debug('Reinstating route_pool database pointers... ')
		# self.route_pool.set_db_pointers(self.db)

	def load_data(self):
		# logger.reading(self.data_path)
		# self._route_pool = pickle.load(open(self.data_path, 'rb'))
		# logger.debug('Reinstating route_pool database pointers... ')
		# self.route_pool.set_db_pointers(self.db)
		raise NotImplementedError

	def generate(self, 
		budget: float = 10000,
		currency: str = 'EUR', 
		max_products = 1000,
		max_reactions = 1000,
		debug=True, 
		max_iter=150, 
		# pick_inner_cheapest=True, 
		# add_size=1,
		shuffle=True,
	):

		from .price import Price

		budget = Price(budget, currency)

		# logger.warning('Something went wrong')

		# raise NotImplementedError


		# assert pick_inner_cheapest
		# assert add_size == 1
		
		# assert len(self.suppliers) == 1

		# supplier = self.suppliers[0]

		recipe = self.starting_recipe.copy()

		recipe.reactants._supplier = self.suppliers

		pool = self.route_pool.copy()

		# return pool

		# pool = list(chain.from_iterable(pool.values()))

		if shuffle:
			from random import shuffle
			logger.debug('Shuffling Route pool')
			shuffle(pool)

		logger.var('route pool', len(pool))
		logger.var('max_iter', max_iter)

		# randomly add products until budget exceeded
		pbar = tqdm()
		# with tqdm() as pbar:
		# with tqdm(total=max_iter) as pbar:
		
		# for i in tqdm(range(max_iter), total=max_iter):
		for i in range(max_iter):

			# logger.title(f'Iteration {i}')
			
			if debug: logger.title(f'Iteration {i}')
			
			price = recipe.price

			if debug: logger.var('price', price)

			# pop a route
			candidate_route = pool.pop()

			if debug: logger.var('candidate_route', candidate_route)
			if debug: logger.var('candidate_route.reactants', candidate_route.reactants.ids)

			# add the route to the recipe
			if debug: logger.var('#recipe.reactants', len(recipe.reactants))
			recipe += candidate_route
			if debug: logger.var('#recipe.reactants', len(recipe.reactants))
			
			# calculate the new price
			new_price = price

			if debug: logger.var('new price', new_price)

			if not len(pool):
				logger.info('Product pool depleted')
				break

			# check breaking conditions
			if new_price > budget:
				logger.info('Budget exceeded')
				recipe = old_recipe.copy()
				continue

			if len(recipe.reactions) > max_reactions:
				pbar.close()
				logger.info('Max #reactions exceeded')
				# recipe = old_recipe.copy()
				break

			if len(recipe.products) > max_products:
				pbar.close()
				logger.info('Max #products exceeded')
				# recipe = old_recipe.copy()
				break
			
			# accept change
			old_recipe = recipe.copy()

			pbar.update(1)
			pbar.set_postfix(dict(price=str(price)))

		else:
			pbar.close()

		### recalculate the products to see if any extra can be had for free?

		# print(recipe)
		# print(recipe.price)

		# new_recipe = Recipe.from_reactants(reactants=recipe.reactants)

		# # calculate fingerprint?

		# print(new_recipe)
		# print(new_recipe.price)

		logger.success(f'Completed after {i} iterations')

		return recipe

	### DUNDERS
	
	def __call__(self, *args, **kwargs):
		return self.generate(*args, **kwargs)
