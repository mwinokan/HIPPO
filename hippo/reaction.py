
import mcol

from .compound import Compound
from .recipe import Recipe

from mlog import setup_logger
logger = setup_logger('HIPPO')

class Reaction:

	"""A :class:`.Reaction` is a simplified representation of a synthetic pathway to create a product :class:`.Compound`. Reactants (also :class:`.Compound` objects) as well as a reaction type are required.

	:class:`.Reaction` objects should not be created directly. Instead use :meth:`.HIPPO.register_reaction` or :meth:`.HIPPO.reactions`
	"""

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
	def id(self) -> int:
		"""Returns the reaction ID"""
		return self._id

	@property
	def type(self) -> str:
		"""Returns the reaction tyoe"""
		return self._type

	@property
	def product(self):
		"""Returns the reaction's product :class:`.Compound`"""
		if isinstance(self._product, int):
			self._product = self.db.get_compound(id=self._product)
		return self._product

	@property
	def product_yield(self) -> float:
		"""Returns the reaction's product yield (fraction)"""
		return self._product_yield

	@property
	def db(self):
		"""Returns a pointer to the parent database"""
		return self._db

	@property
	def reactants(self):
		"""Returns a :class:`.CompoundSubset` of the reactants"""
		from .cset import CompoundSubset
		return CompoundSubset(self.db, indices=self.reactant_ids)

	@property
	def reaction_str(self) -> str:
		"""Returns a string representing the reaction"""
		s = ' + '.join([str(r) for r in self.reactants])
		s = f'{s} -> {str(self.product)}'
		return s

	@property
	def reactant_ids(self) -> set:
		"""Returns a set of reactant ID's"""
		return set(v for v in self.get_reactant_ids())

	@property
	def product_id(self) -> int:
		"""Returns the product ID"""
		return self.product.id
	
	### METHODS

	def get_reactant_amount_pairs(self) -> list[Compound]:
		"""Returns pairs of reactants and their amounts"""

		compound_ids = self.db.select_where(query='reactant_compound, reactant_amount', table='reactant', key='reaction', value=self.id, multiple=True)

		if compound_ids:
			reactants = [(self.db.get_compound(id=id),amount) for id, amount in compound_ids]
		else:
			return []

		return reactants

	def get_reactant_ids(self) -> list[Compound]:

		compound_ids = self.db.select_where(query='reactant_compound', table='reactant', key='reaction', value=self.id, multiple=True)

		if compound_ids:
			return [id for id, in compound_ids]
		else:
			return []

	def get_ingredients(self, amount, return_reactions=False, return_intermediates=False):
		"""recursively assemble a list of ingredients and reactions required to make the compound"""

		ingredients = []
		intermediates = []
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
				_ingredients, _reactions, _intermediates = reaction.get_ingredients(reactant_amount, return_reactions=True, return_intermediates=True)

				ingredients += _ingredients
				reax += _reactions
				intermediates += _intermediates

				intermediates.append(reaction.product.as_ingredient(reactant_amount))

				reax.append(reaction)

			else:
				ingredient = reactant.as_ingredient(reactant_amount)
				ingredients.append(ingredient)

		output = [ingredients]

		if return_reactions:
			reax.append(self)
			output.append(reax)

		if return_intermediates:
			output.append(intermediates)

		match len(output):
			case 1:
				return output[0]
			case 2:
				return output[0], output[1]
			case 3:
				return output[0], output[1], output[2]

	def get_recipe(self, 
		amount: float = 1, # in mg
	):
		"""Get a :class:`.Recipe` describing how to make the product"""

		products = [self.product.as_ingredient(amount=amount)]

		reactants, reactions, intermediates = self.get_ingredients(amount=amount, return_reactions=True, return_intermediates=True)

		recipe = Recipe(products=products, reactants=reactants, reactions=reactions, intermediates=intermediates)
		
		return recipe

	def summary(self, amount=1):
		"""Print a summary of this reaction's information"""

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

	def draw(self):
		"""Draw the molecules involved in this reaction"""

		from molparse.rdkit import draw_grid

		reactants = self.reactants

		product = self.product

		mols = [r.mol for r in reactants]
		mols.append(product.mol)

		labels = [f'+ {r}' if i>0 else f'{r}' for i,r in enumerate(reactants)]
		labels.append(f'-> {product}')

		return draw_grid(mols, labels=labels, highlightAtomLists=None)

	### DUNDERS

	def __str__(self):
		return f'R{self.id}'

	def __repr__(self):
		return f'{mcol.bold}{mcol.underline}{self}: {self.reaction_str} via {self.type}{mcol.unbold}{mcol.ununderline}'

	def __eq__(self, other):

		if self.type != other.type:
			return False

		if self.product != other.product:
			return False

		if self.reactant_ids != other.reactant_ids:
			return False

		return True
