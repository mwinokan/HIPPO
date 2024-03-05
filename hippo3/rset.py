
from .db import Database

import mcol

import os

import logging
logger = logging.getLogger('HIPPO')

from .reaction import Reaction

class ReactionSet:

	def __init__(self, 
		db: Database, 
		table: str = 'reaction',
	) -> None:
		
		self._db = db
		self._table = table

	### PROPERTIES

	@property
	def db(self):
		return self._db

	@property
	def table(self):
		return self._table

	### DUNDERS

	def __getitem__(self, key) -> Reaction:
		
		match key:

			case int():

				if key == 0:
					return self.__getitem__(key=1)

				if key < 0:
					key = len(self) + 1 + key
					return self.__getitem__(key=key)

				else:
					return self.db.get_reaction(table=self.table, id=key)

			case str():
				return self.db.get_reaction(table=self.table, name=key)

			case _:
				logger.error(f'Unsupported type for ReactionSet.__getitem__(): {key=} {type(key)}')

		return None

	def __repr__(self) -> str:
		return f'{mcol.bold}{mcol.underline}set(R x {len(self)}){mcol.unbold}{mcol.ununderline}'

	def __len__(self) -> int:
		return self.db.count(self.table)

	def __iter__(self):
		return iter(self[i+1] for i in range(len(self)))
