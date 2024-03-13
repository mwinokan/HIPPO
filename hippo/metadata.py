
from collections import UserDict
from typing import Mapping

class MetaData(UserDict):

	def __init__(self, 
		__dict: Mapping[str, str]
	) -> None:

		super().__init__()
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
		self.data.update(data)
		self._update_db(commit=commit)
