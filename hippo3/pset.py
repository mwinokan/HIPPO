
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

	@property
	def names(self):
		result = self.db.select(table=self.table, query='pose_name', multiple=True)
		return [q for q, in result]

	### METHODS

	def get_by_tag(self,tag):
		values = self.db.select_where(query='tag_pose', table='tag', key='name', value=tag, multiple=True)
		ids = [v for v, in values if v]
		return self[ids]

	### DUNDERS

	def __getitem__(self, key) -> Pose:
		
		match key:

			case int():
				if key == 0:
					return self.__getitem__(key=1)

				if key < 0:
					key = len(self) + 1 + key
					return self.__getitem__(key=key)

				else:
					return self.db.get_pose(table=self.table, id=key)

			case str():
				if len(key.split(' ')) > 1:
					return self.db.get_pose(longname=key)
				else:
					return self.db.get_pose(name=key)

			case list():
				return PoseSubset(self.db, self.table, key)

			case tuple():
				return PoseSubset(self.db, self.table, key)

			case slice():

				start = key.start or 1
				stop = key.stop or len(self)
				step = key.step or 1

				indices = [i for i in range(start, stop, step)]

				return self[indices]

			case _:
				logger.error(f'Unsupported type for PoseSet.__getitem__(): {type(key)}')

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

	### PROPERTIES

	@property
	def indices(self):
		return self._indices

	@property
	def names(self):
		return [self.db.select_where(table=self.table, query='pose_name', key='id', value=i, multiple=False)[0] for i in self.indices]

	### DUNDERS

	def __repr__(self) -> str:
		return f'{mcol.bold}{mcol.underline}subset(P x {len(self)}){mcol.unbold}{mcol.ununderline}'

	def __len__(self) -> int:
		return len(self.indices)

	def __iter__(self):
		return iter(self.db.get_pose(table=self.table, id=i) for i in self.indices)

	def __getitem__(self, key) -> Pose:
		try:
			index = self.indices[key]
		except IndexError:
			logger.exception(f'list index out of range: {key=} for {self}')
			raise
		return self.db.get_pose(table=self.table, id=index)
