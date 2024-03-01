
from collections.abc import MutableSet

class TagSet(MutableSet):

	def __init__(self, 
		parent, # Compound or Pose
		tags=(), 
		immutable=False
	):
		
		self._elements = []
		self._immutable = immutable
		self._parent = parent

		for tag in tags:
			self.add(tag)
		
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

	def _add_tag_to_db(self, tag):
		payload = { 'name':tag, self.parent.table:self.parent.id }
		self.db.insert_tag(**payload)
			

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

	def add(self, tag):
		assert not self.immutable
		if tag not in self._elements:
			self._elements.append(tag)
			self._add_tag_to_db(tag)

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
