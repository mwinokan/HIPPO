
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

	### DUNDERS

	def __getitem__(self, key) -> Pose:
		
		match key:

			case int():
				return self.db.get_pose(id=key)

			case str():
				return self.db.get_pose(longname=key)

			case _:
				logger.error(f'Unsupported type for CompoundSet.__getitem__(): {type(key)}')

		return None

	def __repr__(self) -> str:
		# return f'PoseSet(table="{self.table}")'
		return f'{mcol.bold}{mcol.underline}set(P x {len(self)}){mcol.unbold}{mcol.ununderline}'

	def __len__(self) -> int:
		return self.db.count(self.table)