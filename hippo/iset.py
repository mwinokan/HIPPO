
import mcol

import logging
logger = logging.getLogger('HIPPO')

class InteractionSet:

	"""
	"""

	_table = 'interaction'

	def __init__(self,
		db: 'Database',
		indices: list = None,
	):

		self._db = db

		indices = indices or []

		if not isinstance(indices, list):
			indices = list(indices)

		indices = [int(i) for i in indices]

		self._indices = sorted(list(set(indices)))
		self._df = None


	### FACTORIES

	@classmethod
	def from_pose(cls, 
		pose: 'Pose | PoseSet',
	) -> 'InteractionSet':

		self = cls.__new__(cls)

		### get the ID's

		from .pset import PoseSet

		if isinstance(pose, PoseSet):

			# check if all poses have fingerprints
			has_invalid_fps, = pose.db.select_where(query='COUNT(1)', table='pose', key=f'pose_id IN {pose.str_ids} AND pose_fingerprint = 0')

			if has_invalid_fps:
				logger.warning(f'{has_invalid_fps} Poses have not been fingerprinted')

			sql = f"""
			SELECT interaction_id FROM interaction
			WHERE interaction_pose IN {pose.str_ids}
			"""
		
		else:

			sql = f"""
			SELECT interaction_id FROM interaction
			WHERE interaction_pose = {pose.id}
			"""

		ids = pose.db.execute(sql).fetchall()

		ids = [i for i, in ids]

		self.__init__(pose.db, ids)
		
		return self

	@classmethod
	def from_residue(cls, 
		db: 'Database',
		residue_number: int,
		chain: None | str = None,
		target: "Target | int" = 1,
	) -> 'InteractionSet':

		"""Get the set of interactions for a given residue number (and chain)
		
		:param db: HIPPO :class:`.Database`
		:param residue_number: the residue number
		:param chain: the chain name / letter, defaults to any chain
		:param target: the protein :class:`.Target` object or ID, defaults to first target in database
		:returns: a :class:`.InteractionSet` object
		"""

		from .target import Target

		self = cls.__new__(cls)

		if isinstance(target, Target):
			target = target.id

		sql = f"""
		SELECT interaction_id FROM interaction
		INNER JOIN feature
		ON interaction_feature = feature_id
		WHERE feature_target = {target}
		AND feature_residue_number = {residue_number}
		"""

		if chain:
			sql += f' AND feature_chain_name = "{chain}"'

		ids = db.execute(sql).fetchall()

		ids = [i for i, in ids]

		self.__init__(db, ids)
		
		return self

	### PROPERTIES

	@property
	def indices(self) -> list[int]:
		return self._indices

	@property
	def ids(self) -> list[int]:
		return self._indices

	@property
	def db(self) -> 'Database':
		"""The associated HIPPO :class:`.Database`"""
		return self._db

	@property
	def table(self) -> str:
		"""Get the name of the database table"""
		return self._table

	@property
	def str_ids(self) -> str:
		"""Return an SQL formatted tuple string of the :class:`.Interaction` IDs"""
		return str(tuple(self.ids)).replace(',)',')')

	@property
	def classic_fingerprint(self) -> dict:
		"""Classic HIPPO fingerprint dictionary, mapping protein :class:`.Feature` ID's to the number of corresponding ligand features (from any :class:`.Pose`)"""
		return self.get_classic_fingerprint()

	@property
	def df(self):

		if not self._df:
			import json
			from pandas import DataFrame
			from molparse.rdkit.features import INTERACTION_TYPES

			records = self.db.select_all_where(table='interaction', key=f'interaction_id IN {self.str_ids}', multiple=True)

			data = []
			for record in records:

				id, feature_id, pose_id, family, atom_ids, prot_coord, lig_coord, distance, energy = record
				
				feature = self.db.get_feature(id=feature_id)
				
				d = dict(id=id)

				d['feature_id'] = feature_id
				d['pose_id'] = pose_id
				d['target_id'] = feature.target
				
				d['type'] = INTERACTION_TYPES[(feature.family, family)]
				
				d['prot_family'] = feature.family
				d['lig_family'] = family
				
				d['residue_name'] = feature.residue_name
				d['residue_number'] = feature.residue_number
				d['chain_name'] = feature.chain_name

				d['distance'] = distance
				d['energy'] = energy
				
				d['prot_coord'] = json.loads(prot_coord)
				d['lig_coord'] = json.loads(lig_coord)

				d['prot_atoms'] = feature.atom_names
				d['lig_atoms'] = atom_ids
				
				data.append(d)

			df = DataFrame.from_records(data=data)

			self._df = df
		
		return self._df
	


	### METHODS

	def summary(self) -> None:
		"""Print a summary of this :class:`.InteractionSet`"""

		logger.header(self)

		for interaction in self:
			# print(interaction)

			# logger.var(f'{interaction.family_str}', f'{interaction.distance:.1f}')
			logger.var(f'{interaction.type} [{interaction.feature.chain_res_name_number_str}]', f'{interaction.distance:.1f}', dict(unit='Å', append=interaction.feature.chain_res_name_number_str))

	def get_classic_fingerprint(self) -> dict:
		"""Classic HIPPO fingerprint dictionary, mapping protein :class:`.Feature` ID's to the number of corresponding ligand features (from any :class:`.Pose`)"""

		pairs = self.db.execute(f"""
		SELECT interaction_feature, COUNT(1) FROM interaction
		WHERE interaction_id IN {self.str_ids}
		GROUP BY interaction_feature
		""").fetchall()

		return {f:c for f,c in pairs}

	### DUNDERS

	def __len__(self) -> int:
		"""The number of interactions in this set"""
		return len(self.indices)

	def __repr__(self) -> str:
		return f'{mcol.bold}{mcol.underline}''{'f'I x {len(self)}''}'f'{mcol.unbold}{mcol.ununderline}'

	def __iter__(self):
		"""Iterate through interactions in this set"""
		return iter(self.db.get_interaction(id=i) for i in self.indices)

	def __getitem__(self, key) -> 'Interaction | InteractionSet':
		"""Get interaction or subsets thereof from this set"""
		match key:
			case int():
				index = self.indices[key]
				return self.db.get_interaction(id=index)
			
			case slice():
				indices = self.indices[key]
				return InteractionSet(self.db, indices)

			case _:
				raise NotImplementedError	
