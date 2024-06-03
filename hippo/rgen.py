
from .recipe import Recipe
from .cset import CompoundSet, IngredientSet

import logging
logger = logging.getLogger('HIPPO')

class RandomRecipeGenerator:

	def __init__(self, db, *, 
		budget: float = 10000,
		currency: str = 'EUR', 
		max_lead_time = None,
		max_reactions = None,
		suppliers: list | None = None,
		start_with: Recipe | CompoundSet | IngredientSet = None,
	):

		logger.debug('RandomRecipeGenerator.__init__()')

		self._db = db
		self._budget = budget
		self._max_lead_time = max_lead_time
		self._suppliers = suppliers
		self._max_reactions = max_reactions

		logger.var('budget', budget)
		logger.var('max_lead_time', max_lead_time)
		logger.var('suppliers', suppliers)

		match start_with:
			case Recipe():
				logger.debug('Starting with provided Recipe')
				recipe = start_with.copy()
				logger.var('starting price', recipe.price[0], dict(unit=recipe.price[1]))
			case CompoundSet():
				logger.debug('Starting with provided CompoundSet')
				raise NotImplementedError
			case IngredientSet():
				logger.debug('Starting with provided IngredientSet')
				raise NotImplementedError
			case None:
				logger.debug('Starting with empty Recipe')
				raise NotImplementedError
			case _:
				logger.error(f'Unrecognised {type(start_with)=}. Restart kernel?')
				raise TypeError

		self._starting_recipe = recipe

		self._reactant_pool = self.get_reactant_pool()

		logger.var('reactant_pool', self.reactant_pool)

		self._product_pool = self.get_product_pool()
		
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
	
	### METHODS

	def get_reactant_pool(self):

		# get all quoted compounds with acceptable supplier and leadtime

		if self.max_lead_time:

			raise NotImplementedError

		elif self.suppliers:

			sql = f'''
			SELECT DISTINCT quote_compound FROM quote
			WHERE quote_supplier IN {self.suppliers_str} 
			'''

		else:

			sql = f'''
			SELECT DISTINCT quote_compound FROM quote 
			'''

		compound_ids = self.db.execute(sql).fetchall()
		compound_ids = [c for c, in compound_ids]

		return CompoundSet(self.db, compound_ids)

	def get_product_pool(self):

		reactants = self.get_reactant_pool()

		recipe = Recipe.from_reactants(reactants=reactants)
		
		return recipe.products

	def generate(self, debug=True):

		recipe = self.starting_recipe.copy()

		# randomly add products until budget exceeded

		raise NotImplementedError

	### DUNDERS
	
	def __call__(self, *args, **kwargs):
		return self.generate(*args, **kwargs)
