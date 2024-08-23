
import mcol

class Target:
	
	"""Object representing a protein target

	.. attention::
	
		:class:`.Target` objects should not be created directly. Instead use :meth:`.HIPPO.register_target` or :meth:`.Pose.target`.
	
	"""

	def __init__(self, 
		db: 'Database', 
		id: int, 
		name: str,
	) -> None:

		self._db = db
		self._id = id
		self._name = name

	@property
	def db(self) -> 'Database':
		"""Returns a pointer to the parent database"""
		return self._db

	@property
	def id(self) -> int:
		"""Returns the target's ID"""
		return self._id

	@property
	def name(self) -> str:
		"""Returns the target's name"""
		return self._name

	@property
	def feature_ids(self) -> list[int]:
		"""Returns the target's feature ID's"""
		feature_ids = self.db.select_where(
			query='feature_id', 
			table='feature', 
			key='target', 
			value=self.id, 
			none=False, 
			multiple=True,
			sort='feature_chain_name, feature_residue_number'
		)

		if feature_ids:
			feature_ids = [v for v, in feature_ids]

		return feature_ids

	@property
	def features(self) -> list['Feature']:
		"""Returns the target's features"""
		if feature_ids := self.feature_ids:
			return [self.db.get_feature(id=i) for i in feature_ids]
		return None

	def __str__(self) -> str:
		"""Unformatted string representation"""
		return f'T{self.id}'

	def __repr__(self) -> str:
		"""Formatted string representation"""
		return f'{mcol.bold}{mcol.underline}{self} "{self.name}"{mcol.unbold}{mcol.ununderline}'

	def calculate_features(self, 
		protein: 'mp.System'
	) -> list['Feature']:

		"""Calculate features from a protein system

		:param protein: `molparse.System` object, likely from :meth:`.Pose.protein_system`
		:returns: a list of :class:`.Feature` objects

		"""

		features = protein.get_protein_features()

		for f in features:
			self.db.insert_feature(
				family=f.family,
				target=self.id,
				atom_names=[a.name for a in f.atoms],
				residue_name=f.res_name,
				residue_number=f.res_number,
				chain_name=f.res_chain,
				commit=False,
			)

		self.db.commit()

		return self.features
