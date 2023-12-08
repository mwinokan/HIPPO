
# set of Poses

from collections.abc import MutableSet

class PoseSet(MutableSet):

	def __init__(self, poses=(), immutable=False):
		
		self._elements = []
		self._immutable = immutable

		for pose in poses:
			self.add(pose)
		
	### FACTORIES

	### PROPERTIES

	@property
	def immutable(self):
		return self._immutable
	
	@immutable.setter
	def immutable(self,b):
		self._immutable = b

	@property
	def poses(self):
		return self._elements

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

	def add(self, compound):
		assert not self.immutable
		if compound not in self._elements:
			self._elements.append(compound)
		else:
			raise ValueError(f'{compound} already in {self}')

	### DUNDERS

	def __contains__(self, pose):
		if isinstance(pose,str):
			return pose in [c.name for c in self.poses]
		else:
			return pose in self.poses

	def __repr__(self):
		if not self:
			return f'PoseSet(empty)'
		else:
			return f'PoseSet(#poses={len(self)}, [{", ".join(p.name for p in self)}])'

	def __len__(self):
		return len(self._elements)

	def __iter__(self):
		return iter(self._elements)

	def __add__(self, other):
		if isinstance(other, PoseSet):
			return PoseSet(self._elements + other._elements)
		elif isinstance(other, list):
			return PoseSet(self._elements + other)

	def __iadd__(self, other):
		if isinstance(other, PoseSet):
			return PoseSet(self._elements + other._elements)
		elif isinstance(other, list):
			return PoseSet(self._elements + other)

	def __getitem__(self, key):
		return self._elements[key]
		
	def __iter__(self):
		return iter(self.poses)
