
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
		# 		logger.var('starting price', recipe.price[0], dict(unit=recipe.price[1]))
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
		# currency: str = 'EUR', 
		max_products = 1000,
		max_reactions = 1000,
		debug=True, 
		maxiter=150, 
		pick_inner_cheapest=True, 
		add_size=1,
		shuffle=True,
	):

		assert pick_inner_cheapest
		assert add_size == 1
		assert len(self.suppliers) == 1

		supplier = self.suppliers[0]

		recipe = self.starting_recipe.copy()

		if shuffle:
			print('shuffling')
			pool = self.product_pool.shuffled()
		else:
			pool = self.product_pool

		# randomly add products until budget exceeded
		for i in tqdm(range(maxiter), total=maxiter):
			# if debug: logger.title(f'Iteration {i}')
			# if debug: logger.var('price', recipe.price[0], dict(unit=recipe.price[1]))

			c_id = pool.pop_id()
			candidates = CompoundSet(self.db, [c_id])

			if debug: logger.var('candidates', candidates.ids)

			candidate_recipe = candidates.get_recipes(
				pick_cheapest=pick_inner_cheapest, 
				supplier=supplier, 
				unavailable_reaction='quiet',
				reaction_checking_cache=self._reaction_checking_cache,
				reaction_reactant_cache=self._reaction_reactant_cache,
			)

			# if debug: logger.var('candidate_recipe', candidate_recipe)
			# if debug: logger.var('candidate_recipe.reactants', candidate_recipe.reactants.ids)

			recipe += candidate_recipe

			new_price = recipe.price

			if debug: logger.var('new price', new_price[0], dict(unit=new_price[1]))

			if not len(pool):
				logger.info('Product pool depleted')
				break

			# check breaking conditions
			if new_price[0] > budget:
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

		return recipe

	### DUNDERS
	
	def __call__(self, *args, **kwargs):
		return self.generate(*args, **kwargs)
