
import mcol

import logging
logger = logging.getLogger('HIPPO')

class Interaction:

	"""
	Interaction class
	"""

	_table = 'interaction'

	def __init__(self,
		db: 'Database',
		id: int,
		feature_id: int,
		pose_id: int,
		family: str,
		atom_ids: list[int],
		prot_coord: list[float],
		lig_coord: list[float],
		distance: float,
		energy: float | None,
	) -> 'Interaction':

		# from interaction table
		self._id = id
		self._feature_id = feature_id
		self._pose_id = pose_id
		self._family = family
		self._atom_ids = atom_ids
		self._prot_coord = prot_coord
		self._lig_coord = lig_coord
		self._distance = distance
		self._energy = energy

		# placeholders
		self._pose = None
		self._feature = None

		self._db = db
		
	### FACTORIES

	### PROPERTIES

	@property
	def id(self):
		return self._id

	@property
	def table(self):
		return self._table

	@property
	def db(self):
		return self._db

	@property
	def family(self):
		return self._family

	@property
	def pose_id(self):
		return self._pose_id

	@property
	def pose(self):
		if not self._pose:
			self._pose = self.db.get_pose(id=self.pose_id)
		return self._pose
	
	@property
	def feature_id(self):
		return self._feature_id

	@property
	def feature(self):
		if not self._feature:
			self._feature = self.db.get_feature(id=self.feature_id)
		return self._feature

	@property
	def atom_ids(self):
		return self._atom_ids

	@property
	def prot_coord(self):
		return self._prot_coord

	@property
	def lig_coord(self):
		return self._lig_coord

	@property
	def distance(self):
		return self._distance

	@property
	def energy(self):
		return self._energy

	@property
	def family_str(self):
		return f'{repr(self.feature)} ~ {self.family}'

	@property
	def type(self) -> str:
		"""Interaction type string"""
		from molparse.rdkit.features import INTERACTION_TYPES
		return INTERACTION_TYPES[(self.feature.family, self.family)]
	
	### METHODS

	def summary(self):
		
		logger.header(f'Interaction {self.id}')
		
		logger.var('feature', self.feature)
		logger.var('pose', self.pose)
		logger.var('family', self.family)
		logger.var('atom_ids', self.atom_ids)
		logger.var('prot_coord', self.prot_coord)
		logger.var('lig_coord', self.lig_coord)
		logger.var('distance', self.distance)
		logger.var('energy', self.energy)


	### DUNDERS

	def __str__(self) -> str:
		"""Plain string representation"""
		# return f'P{self.pose_id} -> I{self.id} ({self.family}) {self.feature}'
		return f'I{self.id}'

	def __repr__(self):
		"""Formatted string representation"""
		return f'{mcol.bold}{mcol.underline}{str(self)}{mcol.unbold}{mcol.ununderline}'
