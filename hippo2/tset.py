
# set of Tags

# from .tag import Tag
from collections.abc import MutableSet

class TagSet(MutableSet):

	def __init__(self, tags=(), immutable=False):
		
		self._elements = []
		self._immutable = immutable

		for tag in tags:
			if ';' in tag:
				[self.add(t.strip()) for t in tag.split(';')]
			else:
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

	### METHODS

	def pop(self):
		assert not self.immutable
		return self._elements.pop()

	def discard(self, key):
		assert not self.immutable
		if key in self:
			i = self._elements.index(key)
			del self._elements[i]
		else:
			raise ValueError(f'{key} not in {self}')

	def remove(self, key):
		assert not self.immutable
		if key in self:
			i = self._elements.index(key)
			del self._elements[i]
		else:
			raise ValueError(f'{key} not in {self}')

	def add(self, tag, duplicate_error=True):
		assert not self.immutable
		# if isinstance(tag, str):
		# 	tag = Tag(tag)
		if tag not in self._elements:
			self._elements.append(tag)
		elif duplicate_error:
			raise ValueError(f'{tag} already in {self}')

	### DUNDERS

	def __contains__(self, tag):
		return tag in self.tags
		# return tag.name in [t.name for t in self.tags]

	def __repr__(self):
		return str(self._elements)

	def __len__(self):
		return len(self._elements)

	def __iter__(self):
		return iter(self._elements)

	def __iter__(self):
		return iter(self.tags)

	def __add__(self, other):
		return TagSet(self.tags+other)
		