
# from .tools import df_row_to_dict

from .db import Database
from .pose import Pose

import mcol

import os

import logging
logger = logging.getLogger('HIPPO')

class PoseSet:
	"""Object representing the 'pose' table in the :class:`.Database`."""

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
		"""Returns the associated :class:`.Database`"""
		return self._db
	
	@property
	def table(self) -> str:
		return self._table

	@property
	def names(self):
		"""Returns the names of child poses"""
		result = self.db.select(table=self.table, query='pose_name', multiple=True)
		return [q for q, in result]

	@property
	def ids(self):
		"""Returns the IDs of child poses"""
		result = self.db.select(table=self.table, query='pose_id', multiple=True)
		return [q for q, in result]

	@property
	def tags(self):
		"""Returns the set of unique tags present in this pose set"""
		values = self.db.select_where(table='tag', query='DISTINCT tag_name', key='tag_pose IS NOT NULL', multiple=True)
		return set(v for v, in values)

	### METHODS

	def get_by_tag(self,tag):
		"""Get all child poses with a certain tag"""
		values = self.db.select_where(query='tag_pose', table='tag', key='name', value=tag, multiple=True)
		ids = [v for v, in values if v]
		return self[ids]

	def get_by_metadata(self, key: str, value: str | None = None):
		"""Get all child podrd with by their metadata. If no value is passed, then simply containing the key in the metadata dictionary is sufficient"""
		results = self.db.select(query='pose_id, pose_metadata', table='pose', multiple=True)
		if value is None:
			ids = [i for i,d in results if d and f'"{key}":' in d]
		else:
			if isinstance(value, str):
				value = f'"{value}"'
			ids = [i for i,d in results if d and f'"{key}": {value}' in d]
		return self[ids]		

	def summary(self):
		"""Print a summary of this pose set"""
		logger.header('PoseSet()')
		logger.var('#poses', len(self))
		logger.var('tags', self.tags)

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
				return self.db.get_pose(name=key)

			case list():
				return PoseSubset(self.db, key)

			case tuple():
				return PoseSubset(self.db, key)

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
	"""Object representing a subset of the 'pose' table in the :class:`.Database`."""

	def __init__(self,
		db: Database,
		indices: list = None,
		*,
		table: str = 'pose',
	):
		self._db = db
		self._table = table

		indices = indices or []

		self._indices = indices

	### PROPERTIES

	@property
	def indices(self):
		"""Returns the ids of poses in this set"""
		return self._indices

	@property
	def ids(self):
		"""Returns the ids of poses in this set"""
		return self._indices

	@property
	def names(self):
		"""Returns the names of poses in this set"""
		return [self.db.select_where(table=self.table, query='pose_name', key='id', value=i, multiple=False)[0] for i in self.indices]

	@property
	def tags(self):
		"""Returns the set of unique tags present in this pose set"""
		values = self.db.select_where(table='tag', query='DISTINCT tag_name', key=f'tag_pose in {tuple(self.ids)}', multiple=True)
		return set(v for v, in values)

	@property
	def compounds(self):
		"""Get the compounds associated to this set of poses"""
		from .cset import CompoundSubset
		ids = self.db.select_where(table='pose', query='DISTINCT pose_compound', key=f'pose_id in {tuple(self.ids)}', multiple=True)
		ids = [v for v, in ids]
		return CompoundSubset(self.db, ids)

	@property
	def num_compounds(self):
		"""Count the compounds associated to this set of poses"""
		return len(self.compounds)

	### METHODS

	def get_by_tag(self,tag):
		"""Get all child poses with a certain tag"""
		values = self.db.select_where(query='tag_pose', table='tag', key='name', value=tag, multiple=True)
		ids = [v for v, in values if v and v in self.ids]
		return PoseSubset(self.db, ids)

	def get_by_metadata(self, key: str, value: str | None = None):
		"""Get all child poses with by their metadata. If no value is passed, then simply containing the key in the metadata dictionary is sufficient"""
		results = self.db.select(query='pose_id, pose_metadata', table='pose', multiple=True)
		if value is None:
			ids = [i for i,d in results if d and f'"{key}":' in d and i in self.ids]
		else:
			if isinstance(value, str):
				value = f'"{value}"'
			ids = [i for i,d in results if d and f'"{key}": {value}' in d and i in self.ids]
		return PoseSubset(self.db, ids)		

	def summary(self):
		"""Print a summary of this pose set"""
		logger.header('PoseSubset()')
		logger.var('#poses', len(self))
		logger.var('#compounds', self.num_compounds)
		logger.var('tags', self.tags)


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
