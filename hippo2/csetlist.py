
import mout
from collections import UserList

class CompoundSetList(UserList):

	def __init__(self, inherited_list=list()):
		super(CompoundSetList, self).__init__(inherited_list)

	def __getitem__(self,key):

		if isinstance(key,slice):
			return CompoundSetList(self.data[key])

		if isinstance(key,int):
			return self.data[key]

		if isinstance(key,str):
			for cs in self.data:
				if cs.name == key:
					return cs
			else:
				mout.error(f'No CompoundSet named {key}')

	@property
	def names(self):
		return [cs.name for cs in self.data]

	def append(self,item):
		if item.name in self.names:
			mout.error(f'CompoundSetList already contains a CompoundSet named {item.name}')
			return

		self.data.append(item)
