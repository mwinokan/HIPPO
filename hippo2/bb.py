
from .compound import Compound

class BuildingBlock(Compound):

	def __init__(self, smiles, tags=None):
		super().__init__(smiles, smiles, tags)

		self._amount = None # mg
	
	### PROPERTIES	

	@property
	def amount(self):
		return self._amount

	@amount.setter
	def amount(self, a):
		mout.debug(f'{self}.amount = {a}')
		self._amount = a
