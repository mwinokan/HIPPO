
class Tag:

	def __init__(self, name):
		
		self._name = name
		
	### FACTORIES

	### PROPERTIES

	@property
	def name(self):
		return self._name

	### METHODS

	### DUNDERS

	def __str__(self):
		return self.name

	def __repr__(self):
		return f'Tag("{self.name}")'
