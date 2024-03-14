
# from .db import Database
from collections.abc import MutableSet

class TagSet:
	"""Object representing the 'tag' table in the :class:`.Database`."""

	def __init__(self, 
		db, 
		table: str = 'tag',
	) -> None:
		
		self._db = db
		self._table = table

	### FACTORIES

	### PROPERTIES

	@property
	def db(self):
		"""Returns the associated :class:`.database`"""
		return self._db
	
	@property
	def table(self) -> str:
		return self._table

	@property
	def unique(self) -> set:
		"""Returns a set of tag names contained in the table"""
		values = self.db.select(table=self.table, query='DISTINCT tag_name', multiple=True)
		return set(v for v, in values)

	### METHODS

	### DUNDERS

	def __getitem__(self, key):
		raise NotImplementedError



class TagSubset(MutableSet):
	"""Object representing a subset of the 'tag' table in the :class:`.Database` belonging to a certain :class:`.Compound` or :class:`.Pose`."""

	def __init__(self, 
		parent, # Compound or Pose
		tags=(), 
		immutable=False,
		commit=True,
	):
		
		self._elements = []
		self._immutable = immutable
		self._parent = parent

		for tag in tags:
			self.add(tag, commit=commit)
		
	### FACTORIES

	### PROPERTIES

	@property
	def tags(self) -> list:
		"""Returns the elements in this set"""
		return self._elements

	@property
	def immutable(self):
		"""Returns whether this set is immutable"""
		return self._immutable
	
	@immutable.setter
	def immutable(self,b):
		self._immutable = b

	@property
	def parent(self):
		"""Returns this set of tags parent :class:`.Compound` or :class:`.Pose`."""
		return self._parent

	@property
	def db(self):
		"""Returns the associated :class:`.database`"""
		return self.parent.db
	

	
	### DATABASE

	def _remove_tag_from_db(self, tag):
		sql = f'DELETE FROM tag WHERE tag_name="{tag}" AND tag_{self.parent.table}={self.parent.id}'
		self.db.execute(sql)

	def _add_tag_to_db(self, tag, commit=True):
		payload = { 'name':tag, self.parent.table:self.parent.id }
		self.db.insert_tag(**payload, commit=commit)
			

	### METHODS

	def pop(self):
		"""Pop the last element"""
		assert not self.immutable
		return self._elements.pop()

	def discard(self, tag):
		"""Discard an element"""
		return self.discard(tag)

	def remove(self, tag):
		"""Remove an element"""
		assert not self.immutable
		if tag in self:
			i = self._elements.index(tag)
			del self._elements[i]
			self._remove_tag_from_db(tag)
		else:
			raise ValueError(f'{tag} not in {self}')

	def add(self, tag, commit=True):
		"""Add an element"""
		assert not self.immutable
		if tag not in self._elements:
			self._elements.append(tag)
			self._add_tag_to_db(tag, commit=commit)

	### DUNDERS

	def __contains__(self, tag):
		return tag in self.tags

	def __repr__(self):
		return str(self._elements)

	def __len__(self):
		return len(self._elements)

	def __iter__(self):
		return iter(self._elements)

	def __iter__(self):
		return iter(self.tags)

	def __add__(self, other):
		raise NotImplementedError
		# return TagSet(self.tags+other)
