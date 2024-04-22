
from .db import Database

import mcol

import os

import logging
logger = logging.getLogger('HIPPO')

from .reaction import Reaction

class ReactionTable:
	"""Object representing the 'reaction' table in the :class:`.Database`."""

	def __init__(self, 
		db: Database, 
		table: str = 'reaction',
	) -> None:
		
		self._db = db
		self._table = table

	### PROPERTIES

	@property
	def db(self):
		"""Returns the associated :class:`.Database`"""
		return self._db

	@property
	def table(self):
		return self._table

	@property
	def types(self):
		result = self.db.select(table=self.table, query='reaction_type', multiple=True)
		return [q for q, in result]

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

class ReactionSet(ReactionTable):
	"""Object representing a subset of the 'reaction' table in the :class:`.Database`."""
	
	def __init__(self,
		db: Database,
		indices: list = None,
		*,
		table: str = 'reaction',
	):
		
		self._db = db
		self._table = table
		indices = indices or []

		if not isinstance(indices, list):
			indices = list(indices)

		assert all(isinstance(i, int) for i in indices)

		self._indices = sorted(list(set(indices)))

	### PROPERTIES

	@property
	def indices(self):
		return self._indices

	@property
	def ids(self):
		return self._indices

	### METHODS

	def add(self, r):
		assert r._table == 'reaction'
		if (id := r.id) in self._indices:
			self._indices.append(id)
	
	### DUNDERS

	def __repr__(self) -> str:
		# return f'{mcol.bold}{mcol.underline}subset(R x {len(self)}){mcol.unbold}{mcol.ununderline}'
		return f'{mcol.bold}{mcol.underline}''{'f'R x {len(self)}''}'f'{mcol.unbold}{mcol.ununderline}'

	def __len__(self) -> int:
		return len(self.indices)

	def __iter__(self):
		return iter(self.db.get_reaction(table=self.table, id=i) for i in self.indices)

	def __getitem__(self, key) -> Reaction:
		try:
			index = self.indices[key]
		except IndexError:
			logger.exception(f'list index out of range: {key=} for {self}')
			raise
		return self.db.get_reaction(table=self.table, id=index)

	def __add__(self, other):
		for reaction in other:
			self.add(reaction)
		return self
		
