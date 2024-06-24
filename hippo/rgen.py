
from .recipe import Recipe
from .cset import CompoundSet, IngredientSet

from tqdm import tqdm

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

		self._db = db
		# self._budget = budget
		self._max_lead_time = max_lead_time
		self._suppliers = suppliers
		self._max_reactions = max_reactions

		# caches
		self._reaction_checking_cache = {}
		self._reaction_reactant_cache = {}

		# logger.var('budget', budget)
		logger.var('max_lead_time', max_lead_time)
		logger.var('suppliers', suppliers)

		# match start_with:
		# 	case Recipe():
		# 		logger.debug('Starting with provided Recipe')
		# 		recipe = start_with.copy()
		# 		logger.var('starting price', recipe.price[0])
		# 	case CompoundSet():
		# 		logger.debug('Starting with provided CompoundSet')
		# 		raise NotImplementedError
		# 	case IngredientSet():
		# 		logger.debug('Starting with provided IngredientSet')
		# 		raise NotImplementedError
		# 	case None:
		# 		logger.debug('Starting with empty Recipe')
		# 		raise NotImplementedError
		# 	case _:
		# 		logger.error(f'Unrecognised {type(start_with)=}. Restart kernel?')
		# 		raise TypeError

		self._starting_recipe = start_with

		self._reactant_pool = self.get_reactant_pool()

		logger.var('reactant_pool', self.reactant_pool)

		logger.debug('Solving product pool...')
		self._product_pool = self.get_product_pool()

		if start_with:
			self._product_pool -= start_with.products.compounds
		
		logger.var('product_pool', self.product_pool)

		logger.debug('Solving route pool...')
		self._route_pool = self.get_route_pool()

	### FACTORIES

	### PROPERTIES

	@property
	def starting_recipe(self):
		return self._starting_recipe

	@property
	def db(self):
		return self._db

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
		
	
	### METHODS

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
		products = Recipe.from_reactants(reactants=reactants, debug=False)
		return products

	def get_route_pool(self):

		# logger.debug("Fetching route ID's from database...")
		route_ids = self.db.select_where(table='route', query='route_id', key=f'route_product IN {self.product_pool.str_ids}', multiple=True)

		# routes = []
		# for route_id, in tqdm(route_ids):
		# 	route = self.db.get_route(id=route_id)
		# 	routes.append(route)

		routes = [self.db.get_route(id=route_id) for route_id, in tqdm(route_ids)]

		from .recipe import RouteSet
		return RouteSet(self.db, routes)

		# route_pool = {}
		# for route_id, route_product in pairs:
		# 	if route_product not in route_pool:
		# 		route_pool[route_product] = []
		# 	# else:

		# 	route_pool[route_product].append(route_id)

		# logger.debug("Constructing Route objects...")
		# for route_product_id, route_ids in tqdm(route_pool.items()):

		# 	recipes = []
			
		# 	# if len(route_ids) > 1:
		# 		# logger.warning(f'{route_product_id} has multiple routes!')

		# 	for route_id in route_ids:

		# 		recipe = self.db.get_route(id=route_id)

		# 		recipes.append(recipe)

		# 	route_pool[route_product_id] = recipes

		# if flatten:
		# 	from itertools import chain
		# 	route_pool = list(chain.from_iterable(route_pool.values()))

		# return route_pool

	# def generate_subrecipes(self):

	# 	from .cset import CompoundSet
	# 	import json

	# 	for compound in self.product_pool:

	# 		recipes = CompoundSet(self.db, [compound.id]).get_recipes(
	# 			pick_cheapest=False, 
	# 			warn_multiple_solutions=False,
	# 		)

	# 		for recipe in recipes:

	# 			data = recipe.get_dict(
	# 				price=False,
	# 				reactant_supplier=False,
	# 				database=False,
	# 				timestamp=False,
	# 				compound_ids_only=True,
	# 				products=False,
	# 			)

	# 			return data

	# 			str_payload = json.dumps(data)

	# 			self.db.insert_recipe(product_id=compound.id, str_payload=str_payload, commit=False)

	# 		break

	# 	raise NotImplementedError

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

		# randomly add products until budget exceeded
		with tqdm(total=max_iter) as pbar:
			
			# for i in tqdm(range(max_iter), total=max_iter):
			for i in range(max_iter):

				if debug: logger.title(f'Iteration {i}')
				price = recipe.price

				pbar.update(i)
				pbar.set_postfix(dict(price=str(price)))

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
					logger.info('Budget exceeded, choosing new candidates')
					recipe = old_recipe.copy()
					continue

				if len(recipe.reactions) > max_reactions:
					logger.info('Max #reactions exceeded, choosing new candidates')
					recipe = old_recipe.copy()
					continue

				if len(recipe.products) > max_products:
					logger.info('Max #products exceeded, choosing new candidates')
					recipe = old_recipe.copy()
					continue
				
				# accept change
				old_recipe = recipe.copy()

				pbar.update(i)
				pbar.set_postfix(dict(price=str(price)))

		return recipe

	### DUNDERS
	
	def __call__(self, *args, **kwargs):
		return self.generate(*args, **kwargs)
