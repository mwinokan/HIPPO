
import mout
import mcol
from .compound import Compound
from .cset import CompoundSet

class Reaction:

	def __init__(self, reaction_type, product, reactants, product_amount=None, amounts=None):
		
		self._type = reaction_type
		self._product = product
		self._product_amount = product_amount
		self._reactants = CompoundSet(f"{self.type} reactants",reactants)

		if amounts is not None:
			if isinstance(amounts, list):
				for comp, amount in zip(self.reactants, amounts):
					comp.amount = amount
			else:
				for comp in self.reactants:
					comp.amount = amounts
		
	### FACTORIES

	# @classmethod
	#     def single_step(cls, reaction_type, smiles1, smiles2, amount1=1, amount2=1):
	
	#         self = cls.__new__(cls)
	
	#         self.__init__(...)
			
	#         return self

	### PROPERTIES

	@property
	def reactants(self):
		return self._reactants

	@property
	def num_reactants(self):
		return len(self.reactants)
	
	@property
	def type(self):
		return self._type

	@property
	def product(self):
		return self._product

	@property
	def product_amount(self):
		return self._product_amount

	@property
	def reactant_names(self):
		return [c.name for c in self.reactants]

	@property
	def reactant_smiles(self):
		return [c.smiles for c in self.reactants]

	# @property
	# def building_blocks(self):
	# 	return self._reactants

	@property
	def dict(self):
		return dict(
			type=self.type,
			product_name=self.product.name,
			product_smiles=self.product.smiles,
			num_reactants=self.num_reactants,
			reactant_names=self.reactant_names,
			reactant_smiles=self.reactant_smiles,
		)
	
	
	### METHODS

	def summary(self):

		mout.header(str(self))

		# mout.var('type',self.type)

		mout.out(f'\n{mcol.func}product')
		mout.var(self.product.name, self.product_amount, unit='mg', symbol='x')

		mout.out(f'\n{mcol.func}reactants [#={len(self.reactants)}]')
		for reactant in self.reactants:
			mout.var(reactant.name, reactant.amount, unit='mg', symbol='x')

	### DUNDERS

	def __repr__(self):
		return f'Reaction(type={self.type}, #reactants={self.num_reactants}, product={self.product.name})'

	def __eq__(self, other):

		if self.num_reactants != other.num_reactants:
			return False

		for reactant in self.reactants:
			if reactant not in other.reactants:
				return False

		return True
