
# from .tools import df_row_to_dict

from .compound import Compound
from .db import Database

import mcol

import os

import logging
logger = logging.getLogger('HIPPO')

class CompoundSet:

	# class to access entries in database tables containing compounds

	def __init__(self, 
		db: Database, 
		table: str = 'compound',
	) -> None:
		
		self._db = db
		self._table = table
		
	### FACTORIES

	### PROPERTIES

	@property
	def db(self) -> Database:
		return self._db
	
	@property
	def table(self) -> str:
		return self._table

	### METHODS

	### DUNDERS

	def __getitem__(self, key) -> Compound:
		
		match key:

			case int():
				return self.db.get_compound(table=self.table, id=key)

			case str():
				return self.db.get_compound(table=self.table, name=key)

			case _:
				logger.error(f'Unsupported type for CompoundSet.__getitem__(): {type(key)}')

		return None

	def __repr__(self) -> str:
		# return f'CompoundSet(table="{self.table}"")'
		return f'{mcol.bold}{mcol.underline}set(C x {len(self)}){mcol.unbold}{mcol.ununderline}'


	def __len__(self) -> int:
		return self.db.count(self.table)
