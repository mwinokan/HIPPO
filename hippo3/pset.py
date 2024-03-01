
# from .tools import df_row_to_dict

from .db import Database
from .pose import Pose

import mcol

import os

import logging
logger = logging.getLogger('HIPPO')

class PoseSet:

	def __init__(self,
		db: Database,
		table: str = 'pose',
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

	def get_by_tag(self,tag):
		values = self.db.select_where(query='tag_pose', table='tag', key='name', value=tag, multiple=True)
		ids = [v for v, in values if v]
		return self[ids]

	### DUNDERS

	def __getitem__(self, key) -> Pose:
		
		match key:

			case int():
				return self.db.get_pose(id=key)

			case str():
				return self.db.get_pose(longname=key)

			case list():
				return PoseSubset(self.db, self.table, key)

			case slice():

				start = key.start or 1
				stop = key.stop or len(self)
				step = key.step or 1

				indices = [i for i in range(start, stop, step)]

				return self[indices]

			case _:
				logger.error(f'Unsupported type for CompoundSet.__getitem__(): {type(key)}')

		return None

	def __repr__(self) -> str:
		# return f'PoseSet(table="{self.table}")'
		return f'{mcol.bold}{mcol.underline}set(P x {len(self)}){mcol.unbold}{mcol.ununderline}'

	def __len__(self) -> int:
		return self.db.count(self.table)

	def __iter__(self):
		return iter(self[i+1] for i in range(len(self)))


class PoseSubset(PoseSet):

	def __init__(self,
		db: Database,
		table: str = 'pose',
		indices: list = None,
	):

		self._db = db
		self._table = table

		indices = indices or []

		self._indices = indices

	def __getitem__(self, key) -> Pose:
		raise NotImplementedError

	@property
	def indices(self):
		return self._indices

	def __repr__(self) -> str:
		return f'{mcol.bold}{mcol.underline}subset(P x {len(self)}){mcol.unbold}{mcol.ununderline}'

	def __len__(self) -> int:
		return len(self.indices)

	def __iter__(self):
		return iter(self.db.get_compound(table=self.table, id=i) for i in self.indices)
