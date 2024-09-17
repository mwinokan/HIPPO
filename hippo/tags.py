# from .db import Database
from collections.abc import MutableSet
import mcol


class TagTable:
    """Object representing the 'tag' table in the :class:`.Database`.

    .. attention::

            :class:`.TagTable` objects should not be created directly. Instead use the :meth:`.HIPPO.tags` property.

    """

    _table = "tag"

    def __init__(
        self,
        db: "Database",
    ) -> None:

        self._db = db

    ### FACTORIES

    ### PROPERTIES

    @property
    def db(self) -> "Database":
        """Returns a pointer to the parent database"""
        return self._db

    @property
    def table(self) -> str:
        """Returns the name of the :class:`.Database` table"""
        return self._table

    @property
    def unique(self) -> set[str]:
        """Returns a set of unique tag names contained in the table"""
        values = self.db.select(
            table=self.table, query="DISTINCT tag_name", multiple=True
        )
        return set(v for v, in values)

    ### METHODS

    def delete(self, tag: str) -> None:
        """Delete all assignments for the given tag"""
        self.db.delete_where(table="tag", key="name", value=tag)

    ### DUNDERS

    def __repr__(self) -> str:
        """Formatted representation of this object"""
        return f"{mcol.bold}{mcol.underline}Tags {self.unique}{mcol.clear}"


class TagSet(MutableSet):
    """Object representing a subset of the 'tag' table in the :class:`.Database` belonging to a certain :class:`.Compound` or :class:`.Pose`.

    .. attention::

            :class:`.TagSet` objects should not be created directly. Instead use the :meth:`.Compound.tags` or :meth:`.Pose.tags` property.

    """

    def __init__(
        self,
        parent: "Compound | Pose",
        tags: list | tuple | None = None,
        immutable: bool = False,
        commit: bool = True,
    ):

        self._elements = []
        self._immutable = immutable
        self._parent = parent

        tags = tags or ()

        for tag in tags:
            self.add(tag, commit=False)

        if commit:
            self.db.commit()

    ### FACTORIES

    ### PROPERTIES

    @property
    def tags(self) -> list:
        """Returns the elements in this set"""
        return self._elements

    @property
    def immutable(self) -> bool:
        """Is this set is immutable?"""
        return self._immutable

    @immutable.setter
    def immutable(
        self,
        b: bool,
    ) -> None:
        self._immutable = b

    @property
    def parent(self):
        """Returns this set of tags parent :class:`.Compound` or :class:`.Pose`."""
        return self._parent

    @property
    def db(self) -> "Database":
        """Returns a pointer to the parent database"""
        return self.parent.db

    ### DATABASE

    def _remove_tag_from_db(
        self,
        tag: str,
    ) -> None:
        """Delete a specific tag assignment for the parent :class:`.Compound`/:class:`.Pose`

        :param tag: tag to delete

        """
        sql = f'DELETE FROM tag WHERE tag_name="{tag}" AND tag_{self.parent.table} = {self.parent.id}'
        self.db.execute(sql)

    def _clear_tags_from_db(
        self,
        tag: str,
    ) -> None:
        """Delete all tag assignments for the parent :class:`.Compound`/:class:`.Pose`

        :param tag: tag to delete

        """
        sql = f"DELETE FROM tag WHERE tag_{self.parent.table} = {self.parent.id}"
        self.db.execute(sql)

    def _add_tag_to_db(
        self,
        tag: str,
        commit: bool = True,
    ) -> None:
        """Assign a given tag to the parent

        :param tag: tag to add
        :param commit: commit the changes? (Default value = True)

        """
        payload = {"name": tag, self.parent.table: self.parent.id}
        self.db.insert_tag(**payload, commit=commit)

    ### METHODS

    def pop(self) -> str:
        """Pop the last element"""
        assert not self.immutable
        return self._elements.pop()

    def discard(
        self,
        tag: str,
    ) -> None:
        """Discard an element

        :param tag: tag to discard

        """

        self.discard(tag)

    def clear(self):
        """Clear all tags"""
        self._elements = []
        self._clear_tags_from_db(self)

    def remove(self, tag: str) -> None:
        """Remove an element

        :param tag: tag to remove
        :raises ValueError: if tag is not in set

        """

        assert not self.immutable
        if tag in self:
            i = self._elements.index(tag)
            del self._elements[i]
            self._remove_tag_from_db(tag)
        else:
            raise ValueError(f"{tag} not in {self}")

    def add(
        self,
        tag: str,
        commit: bool = True,
    ) -> None:
        """Add a tag to the set

        :param tag: tag to add
        :param commit: commit the change? (Default value = True)

        """

        assert not self.immutable
        if tag not in self._elements:
            self._elements.append(tag)
            self._add_tag_to_db(tag, commit=commit)

    def glob(self, pattern: str) -> list[str]:
        """Construct a list from tags in the set names that match a given UNIX-style pattern.

        :param pattern: unix style pattern with shell-style wildcards
        :returns: list of tags

        """

        import fnmatch

        return fnmatch.filter(self.tags, pattern)

    ### DUNDERS

    def __contains__(self, tag: str) -> bool:
        """Is this tag in the set?"""
        return tag in self.tags

    def __repr__(self) -> str:
        """Formatted representation of this set"""
        return str(self._elements)

    def __len__(self) -> int:
        """Number of tags in this set"""
        return len(self._elements)

    def __iter__(self):
        """Iterate through this set"""
        return iter(self._elements)

    def __add__(self, other):
        """
        .. attention::

                Adding sets together is not supported

        """
        raise NotImplementedError

    def __getitem__(self, key: int):
        """Get a specific element in the set by index"""
        return self._elements[key]
