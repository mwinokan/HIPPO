
import mout

class Pose:

	def __init__(self, name, compound, pdb_path, site_index, chain=None):

		# arguments
		self._name = name
		self._compound = compound
		self._pdb_path = pdb_path
		self._site_index = site_index
		self._chain = chain

		# blanks
		self._fingerprint = None
		
	### FACTORIES

	@classmethod
	def from_bound_pdb(cls, compound, pdb_path, metadata, site_index=None, chain=None):

		self = cls.__new__(cls)

		pose_name = metadata['crystal_name']

		if site_index is None:
			site_index = pose_name[-2]

		if chain is None:
			chain = pose_name[-1]

		self.__init__(f'{site_index}{chain}', compound, pdb_path, site_index, chain)
		
		return self

	### PROPERTIES

	@property
	def compound(self):
		return self._compound

	@property
	def name(self):
		return self._name
	
	@property
	def fingerprint(self):

		if self._fingerprint is None:
			self.calculate_fingerprint()

		return self._fingerprint

	### METHODS
	
	def calculate_fingerprint(self):
		fingerprint = None
		return fingerprint

	### DUNDERS

	def __repr__(self):
		return f'Pose("{self.name}", compound={self.compound.name})'
