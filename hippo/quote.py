
import mcol
# from dataclasses import dataclass, field, asdict

from .price import Price

import logging
logger = logging.getLogger('HIPPO')

# from typing import TypeVar, Type

# T = TypeVar('T')

# @dataclass
class Quote:
	
	"""Dataclass representing a quote in the database"""

	_db = None

	def __init__(self,
		db,
		id: int,
		compound: int,
		smiles: str,
		supplier: str,
		catalogue: str,
		entry: str,
		amount: float,
		price: float,
		currency: str,
		purity: float,
		lead_time: int,
		date: str | None = None,
		type: str | None = None,
	):

		price = Price(price, currency)

		self._db = db
		self._id = id
		self._compound = compound
		self._smiles = smiles
		self._supplier = supplier
		self._catalogue = catalogue
		self._entry = entry
		self._amount = amount
		self._price = price
		# self._currency = currency
		self._purity = purity
		self._lead_time = lead_time
		self._date = date
		self._type = type

		from datetime import datetime
		quote_age = (datetime.today() - datetime.strptime(self.date, '%Y-%m-%d')).days
		# if quote_age > 30:
			# logger.warning(f'Quote is {quote_age} days old')
			# logger.warning(self)

	### FACTORIES

	@classmethod
	def combination(cls, required_amount, quotes):

		"""

		* Start with biggest pack
		* Estimate by scaling linearly with unit price
		* Use MCule confirmed stock amount

		"""

		biggest_pack = sorted(quotes, key=lambda x: x.amount)[-1]
		unit_price = biggest_pack.price / biggest_pack.amount
		estimated_price = unit_price * required_amount

		self = cls.__new__(cls)
		self.__init__(
			db=biggest_pack.db,
			id=None,
			compound=biggest_pack.compound,
			smiles=biggest_pack.smiles,
			supplier=biggest_pack.supplier,
			catalogue=biggest_pack.catalogue,
			entry=biggest_pack.entry,
			amount=required_amount,
			price=estimated_price,
			currency=biggest_pack.currency,
			purity=biggest_pack.purity,
			lead_time=biggest_pack.lead_time,
			date=biggest_pack.date,
			type=f'estimate from quote={biggest_pack.id}'
		)
		
		return self

	### PROPERTIES

	@property
	def entry_str(self):
		if self.catalogue:
			return f'{self.supplier}:{self.catalogue}:{self.entry}'
		else:
			return f'{self.supplier}:{self.entry}'

	@property
	def db(self):
		return self._db

	@property
	def id(self):
		return self._id

	@property
	def compound(self):
		return self._compound

	@property
	def smiles(self):
		return self._smiles

	@property
	def supplier(self):
		return self._supplier

	@property
	def catalogue(self):
		return self._catalogue

	@property
	def entry(self):
		return self._entry

	@property
	def amount(self):
		return self._amount

	@property
	def price(self):
		return self._price

	@property
	def currency(self):
		return self.price.currency

	@property
	def purity(self):
		return self._purity

	@property
	def lead_time(self):
		return self._lead_time

	@property
	def date(self):
		return self._date

	@property
	def type(self):
		return self._type

	@property
	def dict(self):
		return dict(
			id=self.id,
			compound=self.compound,
			smiles=self.smiles,
			supplier=self.supplier,
			catalogue=self.catalogue,
			entry=self.entry,
			amount=self.amount,
			price=self.price,
			# currency=self.currency,
			purity=self.purity,
			lead_time=self.lead_time,
			date=self.date,
			type=self.type,
		)

	@property
	def currency_symbol(self):
		return self.price.symbol
	
	### METHODS

	### DUNDERS

	def __repr__(self):

		if self.purity:
			purity = f' @ {self.purity:.0%}'
		else:
			purity = ''

		if self.supplier == 'Stock':
			return f'{mcol.bold}{mcol.underline}C{self.compound} In Stock: {self.amount:}mg{purity}{mcol.unbold}{mcol.ununderline}'

		if self.type:
			return f'{mcol.bold}{mcol.underline}C{self.compound} {self.entry_str} {self.amount:}mg{purity} = {self.price:} ({self.lead_time} days) {self.smiles} [{self.type}]{mcol.unbold}{mcol.ununderline}'
		else:
			return f'{mcol.bold}{mcol.underline}C{self.compound} {self.entry_str} {self.amount:}mg{purity} = {self.price:} ({self.lead_time} days) {self.smiles}{mcol.unbold}{mcol.ununderline}'
