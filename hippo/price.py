
import mcol

import logging
logger = logging.getLogger('HIPPO')

CURRENCIES = {
	'USD':'$', 
	'EUR':'€', 
	'GBP':'£',
}

class Price:

	def __init__(self, amount, currency):

		if currency not in CURRENCIES:
			assert currency is None
			assert not amount
			amount = None
		
		self._amount = amount
		self._currency = currency

		# assert currency in CURRENCIES, f"Unrecognised {currency=}"
		
	### FACTORIES

	@classmethod
	def null(cls):
		self = cls.__new__(cls)
		self.__init__(None, None)
		return self

	### PROPERTIES

	@property
	def symbol(self):
		return CURRENCIES[self.currency]
	
	@property
	def currency(self):
		return self._currency
	
	@property
	def amount(self):
		return self._amount

	@property
	def is_null(self):
		return self.amount is None

	### METHODS

	### DUNDERS

	def __str__(self):
		if self.currency is None:
			return 'Null'

		return f'{self.symbol}{self.amount:.2f} {self.currency}'

	def __repr__(self):
		if self.currency is None:
			return f'{mcol.bold}{mcol.underline}Null Price{mcol.unbold}{mcol.ununderline}'

		return f'{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}'

	def __add__(self, other):

		if self.is_null:
			return other
		
		if other.is_null:
			return self

		if self.currency != other.currency:
			raise NotImplementedError(f'Adding two different currencies: {self.currency} != {other.currency}')
		return Price(self.amount + other.amount, self.currency)

	# def __cmp__(self, other):
		# assert self.currency == other.currency
		# return self.amount - other.amount

	def __lt__(self, other):

		if self.is_null and other.is_null:
			return False

		if self.is_null and not other.is_null:
			return True

		if not self.is_null and other.is_null:
			return False

		assert self.currency == other.currency, f'Comparing different currencies: {self.currency} != {other.currency}'
		return self.amount < other.amount

	def __gt__(self, other):

		if self.is_null and other.is_null:
			return False

		if self.is_null and not other.is_null:
			return False

		if not self.is_null and other.is_null:
			return True

		assert self.currency == other.currency, f'Comparing different currencies: {self.currency} != {other.currency}'
		return self.amount > other.amount

# class ComboPrice:

	# ...

