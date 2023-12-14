
from .compound import Compound
import mout
import math

class BuildingBlock(Compound):

	def __init__(self, smiles, tags=None, amount=None):
		super().__init__(smiles, smiles, tags)

		self._amount = amount # mg

		self._name_is_smiles = True

		self._lead_time = None
		self._price_picker = None

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
			has_price_picker = self.price_picker is not None,
			lead_time = self.lead_time,
			quote_attempted = 'quote_attempted' in self.tags,
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

	@property
	def lead_time(self):
		return self._lead_time
	
	@lead_time.setter
	def lead_time(self, d):
		self._lead_time = d

	@property
	def price_picker(self):
		return self._price_picker
	
	@price_picker.setter
	def price_picker(self, d):
		self._price_picker = d

	### METHODS

	def copy(self):
		bb = BuildingBlock(self.smiles, self.tags, amount=self.amount)
		bb._name = self.name
		bb._name_is_smiles = self.name_is_smiles
		return bb

	def get_price(self, *args, **kwargs):
		assert self.price_picker
		return self.price_picker.get_price(*args, **kwargs)

	def get_pack(self, *args, **kwargs):
		assert self.price_picker
		return self.price_picker.get_pack(*args, **kwargs)

	### DUNDERS

	def __repr__(self):
		return f'BuildingBlock({self.name}, {self.smiles}, #amount={self.amount})'

class PricePicker:

	def __init__(self, price_dict):

		self._data = {}
		for d in sorted(price_dict, key=lambda x: x['amount']):
			self._data[d['amount']] = d['price']

		self._min_amount = min(self.data.keys())
		self._min_price = self.data[self.min_amount]
		self._max_amount = max(self.data.keys())
		self._max_price = self.data[self.max_amount]

	def get_pack(self, amount):

		budget = self.get_price(amount)
		order = amount

		for amount, price in self.data.items():
			if price <= budget and amount >= order:
				order = amount

		return dict(amount=order, price=budget)

	def get_price(self, query):
		for amount, price in self.data.items():
			if amount >= query:
				return price

		else:
			return self.max_price * query/self.max_amount

		# if amount < self.min_amount:
			# return self.min_price

	# def __getitem__(self, key):
	# 	return self._data[key]

	@property
	def min_amount(self):
		return self._min_amount

	@property
	def min_price(self):
		return self._min_price

	@property
	def max_amount(self):
		return self._max_amount

	@property
	def max_price(self):
		return self._max_price

	@property
	def data(self):
		return self._data
	