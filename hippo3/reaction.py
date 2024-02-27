
import mcol

from .compound import Compound
from .recipe import Recipe

from mlog import setup_logger
logger = setup_logger('HIPPO', debug=True)

class Reaction:

	def __init__(self,
		db,
		id: int,
		type: str,
		product: int,
		product_yield: float,
	):
		
		self._db = db
		self._id = id
		self._type = type
		self._product = product
		self._product_yield = product_yield
		
	### FACTORIES

	### PROPERTIES

	@property
	def id(self):
		return self._id

	@property
	def type(self):
		return self._type

	@property
	def product(self):
		return self._product

	@property
	def product_yield(self):
		return self._product_yield

	@property
	def db(self):
		return self._db

	### METHODS

	def get_reactant_amount_pairs(self) -> list[Compound]:

		compound_ids = self.db.select_where(query='reactant_compound, reactant_amount', table='reactant', key='reaction', value=self.id, multiple=True)

		if compound_ids:
			reactants = [(self.db.get_compound(id=id),amount) for id, amount in compound_ids]
		else:
			return []

		return reactants

	def get_ingredients(self, amount, return_reactions=False):

		ingredients = []
		reax = []

		pairs = self.get_reactant_amount_pairs()
		  
		for reactant, reactant_amount in pairs:

			# scale amount
			reactant_amount *= amount #/total_reactant_amount
			reactant_amount /= self.product_yield

			reactions = reactant.get_reactions(none='quiet')

			if reactions:
				assert len(reactions) == 1
				reaction = reactions[0]
				reax.append(reaction)
				ingredients += reaction.get_ingredients(reactant_amount)

			else:
				
				ingredient = reactant.as_ingredient(reactant_amount)

				ingredients.append(ingredient)

		if return_reactions:
			reax.append(self)
			return ingredients, reax
			
		return ingredients

	def get_recipe(self, amount):

		ingredients, reactions = self.get_ingredients(amount=amount, return_reactions=True)

		recipe = Recipe(ingredients, reactions)
		
		return recipe

	def summary(self, amount=1):

		print(f'{self}.id={self.id}')
		print(f'{self}.type={self.type}')
		print(f'{self}.product={self.product}')
		print(f'{self}.product_yield={self.product_yield}')

		reactants = self.get_reactant_amount_pairs()
		print(f'{self}.reactants={reactants}')

		print(f'Ingredients for {amount} mg of product:')
		ingredients = self.get_ingredients(amount=amount)
		print(ingredients)

		return self.get_recipe(amount)

	### DUNDERS

	def __repr__(self):
		# return f'Reaction(#{self.id})'
		return f'{mcol.bold}{mcol.underline}R{self.id}{mcol.unbold}{mcol.ununderline}'

