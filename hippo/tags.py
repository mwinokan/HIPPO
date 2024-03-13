
# from .db import Database
from collections.abc import MutableSet

class TagSet:

	# class to access entries in database tables containing tags

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
		return self._db
	
	@property
	def table(self) -> str:
		return self._table

	@property
	def unique(self):
		values = self.db.select(table=self.table, query='DISTINCT tag_name', multiple=True)
		return set(v for v, in values)

	### METHODS

	### DUNDERS

	def __getitem__(self, key):
		raise NotImplementedError



class TagSubset(MutableSet):

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
	def tags(self):
		return self._elements

	@property
	def immutable(self):
		return self._immutable
	
	@immutable.setter
	def immutable(self,b):
		self._immutable = b

	@property
	def parent(self):
		return self._parent

	@property
	def db(self):
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
		assert not self.immutable
		return self._elements.pop()

	def discard(self, tag):
		return self.discard(tag)

	def remove(self, tag):
		assert not self.immutable
		if tag in self:
			i = self._elements.index(tag)
			del self._elements[i]
			self._remove_tag_from_db(tag)
		else:
			raise ValueError(f'{tag} not in {self}')

	def add(self, tag, commit=True):
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
