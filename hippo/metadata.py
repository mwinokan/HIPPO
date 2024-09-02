from collections import UserDict
from typing import Mapping


class MetaData(UserDict):
    """Metadata dictionary linked to a compound or pose in a HIPPO :class:`.Database`

    .. attention::

            :class:`.Metadata` objects should not be created directly. Instead use the methods :meth:`.Compound.metadata` and :meth:`.Pose.metadata`

    """

    def __init__(
        self,
        __dict: Mapping[str, str] | None,
    ) -> None:

        super().__init__()
        if __dict:
            for key, value in __dict.items():
                super().__setitem__(key, value)

        self._db = None
        self._table: str = None
        self._id: str = None

    ### PROPERTIES

    @property
    def table(self) -> str:
        """name of the associated :class:`.Database` table"""
        return self._table

    @property
    def id(self) -> int:
        """entry ID in the associated :class:`.Database` table"""
        return self._id

    @property
    def db(self) -> "Database":
        """associated :class:`.Database`"""
        return self._db

    ### METHODS

    def _update_db(
        self,
        commit: bool = True,
    ) -> None:
        """Update the associated :class:`.Database` entry

        :param commit: commit the changes (Default value = True)

        """
        self._db.insert_metadata(
            table=self._table, id=self._id, payload=self.data, commit=commit
        )

    def update(
        self,
        data: dict,
        commit: bool = True,
    ) -> None:
        """Wrapper for dict.update()

        :param data: data with which to update the metadata
        :param commit: commit the changes (Default value = True)

        """
        self.data.update(data)
        self._update_db(commit=commit)

    def append(
        self,
        key: str,
        value,
        commit: bool = True,
    ) -> None:
        """Create or append to a list-like value with given key

        :param key: metadata dictionary key to be modified
        :param value: value to be appended to list-like ``metadata[key]``
        :param commit: commit the changes (Default value = True)

        """
        if key not in self:
            self.data[key] = []
        if value not in self.data[key]:
            self.data[key].append(value)
        self._update_db(commit=commit)

    ### DUNDERS

    def __setitem__(
        self,
        key: str,
        item,
    ) -> None:
        """Set the value associated to a specific dictionary key"""

        self.data.__setitem__(key, item)
        self._update_db()

    def __delitem__(self, key: str) -> None:
        """Set a dictionary key"""

        self.data.__delitem__(key)
        self._update_db()
