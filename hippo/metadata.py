
from collections import UserDict
from typing import Mapping

class MetaData(UserDict):

	"""Metadata dictionary linked to a compound or pose in the database"""

	def __init__(self, 
		__dict: Mapping[str, str] | None
	) -> None:

		super().__init__()
		if __dict:
			for key, value in __dict.items():
				super().__setitem__(key, value)

		self._db = None
		self._table: str = None
		self._id: str = None

	def _update_db(self, commit=True):
		self._db.insert_metadata(table=self._table, id=self._id, payload=self.data, commit=commit)

	def __setitem__(self, key: str, item: str) -> None:
		self.data.__setitem__(key, item)
		self._update_db()

	def __delitem__(self, key: str) -> None:
		self.data.__delitem__(key)
		self._update_db()

	def update(self, data, commit=True):
		"""Wrapper for dict.update()"""
		self.data.update(data)
		self._update_db(commit=commit)

	def append(self, key, value):
		"""Create or append to a list-like value with given key"""
		if key not in self:
			self.data[key] = []
		if value not in self.data[key]:
			self.data[key].append(value)
		self._update_db()
