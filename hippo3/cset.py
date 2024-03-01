
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

	@property
	def names(self):
		result = self.db.select(table=self.table, query='compound_name', multiple=True)
		return [q for q, in result]

	### METHODS

	### DUNDERS

	def __getitem__(self, key) -> Compound:
		
		match key:

			case int():
				return self.db.get_compound(table=self.table, id=key)

			case str():
				return self.db.get_compound(table=self.table, name=key)

			case slice():

				start = key.start or 1
				stop = key.stop or len(self)
				step = key.step or 1

				indices = [i for i in range(start, stop, step)]

				return CompoundSubset(self.db, self.table, indices)

			case _:
				logger.error(f'Unsupported type for CompoundSet.__getitem__(): {key=} {type(key)}')

		return None

	def __repr__(self) -> str:
		return f'{mcol.bold}{mcol.underline}set(C x {len(self)}){mcol.unbold}{mcol.ununderline}'

	def __len__(self) -> int:
		return self.db.count(self.table)

	def __iter__(self):
		return iter(self[i+1] for i in range(len(self)))


class CompoundSubset(CompoundSet):

	def __init__(self,
		db: Database,
		table: str = 'compound',
		indices: list = None,
	):

		self._db = db
		self._table = table

		indices = indices or []

		self._indices = indices

		print(indices)

	def __getitem__(self, key) -> Compound:
		raise NotImplementedError

	@property
	def indices(self):
		return self._indices

	def __repr__(self) -> str:
		return f'{mcol.bold}{mcol.underline}subset(C x {len(self)}){mcol.unbold}{mcol.ununderline}'

	def __len__(self) -> int:
		return len(self.indices)

	def __iter__(self):
		return iter(self.db.get_compound(table=self.table, id=i) for i in self.indices)
