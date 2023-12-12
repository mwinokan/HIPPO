
from .compound import Compound
import mout

class BuildingBlock(Compound):

	def __init__(self, smiles, tags=None, amount=None):
		super().__init__(smiles, smiles, tags)

		self._amount = amount # mg

		self._name_is_smiles = True

		# self._enamine_bb_id = None

		## blanks
		# self._purchaseable = None
	
	### PROPERTIES	

	@property
	def amount(self):
		return self._amount

	@amount.setter
	def amount(self, a):
		# mout.debug(f'{self}.amount = {a}')
		self._amount = a
	
	@property
	def dict(self):
		return dict(
			name = self.name,
			smiles = self.smiles,
			amount = self.amount,
			name_is_smiles = self.name_is_smiles,
		)

	@property
	def name(self):
		return self._name

	@name.setter
	def name(self, name):
		self._name = name
		self._name_is_smiles = False

	@property
	def name_is_smiles(self):
		return self._name_is_smiles
	
	### METHODS

	def copy(self):
		bb = BuildingBlock(self.smiles, self.tags, amount=self.amount)
		bb._name = self.name
		bb._name_is_smiles = self.name_is_smiles
		return bb

	### DUNDERS

	def __repr__(self):
		return f'BuildingBlock({self.name}, {self.smiles}, #amount={self.amount})'
