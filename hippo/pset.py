
# from .tools import df_row_to_dict

from .db import Database
from .pose import Pose

import mcol

import os

import logging
logger = logging.getLogger('HIPPO')

class PoseTable:
	"""Object representing the 'pose' table in the :class:`.Database`."""

	def __init__(self,
		db: Database,
		table: str = 'pose',
	) -> None:

		self._db = db
		self._table = table
		
	### FACTORIES

	### PROPERTIES
	
	@property
	def db(self) -> Database:
		"""Returns the associated :class:`.Database`"""
		return self._db
	
	@property
	def table(self) -> str:
		return self._table

	@property
	def names(self):
		"""Returns the names of child poses"""
		result = self.db.select(table=self.table, query='pose_name', multiple=True)
		return [q for q, in result]

	@property
	def ids(self):
		"""Returns the IDs of child poses"""
		result = self.db.select(table=self.table, query='pose_id', multiple=True)
		return [q for q, in result]

	@property
	def tags(self):
		"""Returns the set of unique tags present in this pose set"""
		values = self.db.select_where(table='tag', query='DISTINCT tag_name', key='tag_pose IS NOT NULL', multiple=True)
		return set(v for v, in values)

	### METHODS

	def get_by_tag(self, tag):
		"""Get all child poses with a certain tag"""
		values = self.db.select_where(query='tag_pose', table='tag', key='name', value=tag, multiple=True)
		ids = [v for v, in values if v]
		# print(values)
		return self[ids]

	def get_by_metadata(self, key: str, value: str | None = None):
		"""Get all child podrd with by their metadata. If no value is passed, then simply containing the key in the metadata dictionary is sufficient"""
		results = self.db.select(query='pose_id, pose_metadata', table='pose', multiple=True)
		if value is None:
			ids = [i for i,d in results if d and f'"{key}":' in d]
		else:
			if isinstance(value, str):
				value = f'"{value}"'
			ids = [i for i,d in results if d and f'"{key}": {value}' in d]
		return self[ids]			

	def summary(self):
		"""Print a summary of this pose set"""
		logger.header('PoseTable()')
		logger.var('#poses', len(self))
		logger.var('tags', self.tags)

	### DUNDERS

	def __getitem__(self, key) -> Pose:
		
		match key:

			case int():
				if key == 0:
					return self.__getitem__(key=1)

				if key < 0:
					key = len(self) + 1 + key
					return self.__getitem__(key=key)

				else:
					return self.db.get_pose(table=self.table, id=key)

			case str():
				return self.db.get_pose(name=key)

			case key if isinstance(key, list) or isinstance(key, tuple) or isinstance(key, set):

				indices = []
				for i in key:
					if isinstance(i,int):
						index = i
					elif isinstance(i,str):
						index = self.db.get_pose_id(name=i)
					else:
						raise NotImplementedError

					assert index
					indices.append(index)

				return PoseSet(self.db, indices)

			case slice():

				start = key.start or 1
				stop = key.stop or len(self)
				step = key.step or 1

				indices = [i for i in range(start, stop, step)]

				return self[indices]

			case _:
				logger.error(f'Unsupported type for PoseTable.__getitem__(): {type(key)}')

		return None

	def __repr__(self) -> str:
		# return f'PoseTable(table="{self.table}")'
		return f'{mcol.bold}{mcol.underline}set(P x {len(self)}){mcol.unbold}{mcol.ununderline}'

	def __len__(self) -> int:
		return self.db.count(self.table)

	def __iter__(self):
		return iter(self[i+1] for i in range(len(self)))


class PoseSet(PoseTable):
	"""Object representing a subset of the 'pose' table in the :class:`.Database`."""

	def __init__(self,
		db: Database,
		indices: list = None,
		*,
		table: str = 'pose',
	):
		self._db = db
		self._table = table

		indices = indices or []

		if not isinstance(indices, list):
			indices = list(indices)

		self._indices = indices

	### PROPERTIES

	@property
	def indices(self):
		"""Returns the ids of poses in this set"""
		return self._indices

	@property
	def ids(self):
		"""Returns the ids of poses in this set"""
		return self._indices

	@property
	def names(self):
		"""Returns the names of poses in this set"""
		return [self.db.select_where(table=self.table, query='pose_name', key='id', value=i, multiple=False)[0] for i in self.indices]

	@property
	def smiles(self):
		"""Returns the smiles of poses in this set"""
		return [self.db.select_where(table=self.table, query='pose_smiles', key='id', value=i, multiple=False)[0] for i in self.indices]

	@property
	def tags(self):
		"""Returns the set of unique tags present in this pose set"""
		values = self.db.select_where(table='tag', query='DISTINCT tag_name', key=f'tag_pose in {tuple(self.ids)}', multiple=True)
		return set(v for v, in values)

	@property
	def compounds(self):
		"""Get the compounds associated to this set of poses"""
		from .cset import CompoundSet
		ids = self.db.select_where(table='pose', query='DISTINCT pose_compound', key=f'pose_id in {tuple(self.ids)}', multiple=True)
		ids = [v for v, in ids]
		return CompoundSet(self.db, ids)

	@property
	def num_compounds(self):
		"""Count the compounds associated to this set of poses"""
		return len(self.compounds)

	@property
	def df(self):
		"""Get a DataFrame of the poses in this set"""
		return self.get_df(mol=True)

	### FILTERING

	def get_by_tag(self, tag, inverse=False):
		"""Get all child poses with a certain tag"""
		values = self.db.select_where(query='tag_pose', table='tag', key='name', value=tag, multiple=True)
		if inverse:
			matches = [v for v, in values if v]
			ids = [i for i in self.ids if i not in matches]
		else:
			ids = [v for v, in values if v and v in self.ids]
		return PoseSet(self.db, ids)

	def get_by_metadata(self, key: str, value: str | None = None):
		"""Get all child poses with by their metadata. If no value is passed, then simply containing the key in the metadata dictionary is sufficient"""
		results = self.db.select(query='pose_id, pose_metadata', table='pose', multiple=True)
		if value is None:
			ids = [i for i,d in results if d and f'"{key}":' in d and i in self.ids]
		else:
			if isinstance(value, str):
				value = f'"{value}"'
			ids = [i for i,d in results if d and f'"{key}": {value}' in d and i in self.ids]
		return PoseSet(self.db, ids)		

	def get_by_inspiration(self, inspiration: int | Pose, inverse=False):
		"""Get all child poses with with this inspiration."""

		ids = set()

		for pose in self:
			if not inverse:
				for pose_inspiration in pose.inspirations:
					if pose_inspiration == inspiration:
						ids.add(pose.id)
						break

			elif inverse:
				for pose_inspiration in pose.inspirations:
					if pose_inspiration == inspiration:
						break
				else:
					ids.add(pose.id)

		return PoseSet(self.db, ids)

	def get_df(self, skip_no_mol=True, **kwargs):
		"""Get a DataFrame of the poses in this set"""
		
		from pandas import DataFrame

		data = []

		for pose in self:
			d = pose.get_dict(**kwargs)

			if not d['mol']:
				logger.warning(f'Skipping pose with no mol: {d["id"]} {d["name"]}')
				continue
			data.append(d)

		return DataFrame(data)

	### TAGGING

	

	### EXPORTING

	def write_sdf(self, out_path, name_col='name'):

		"""Write an SDF"""

		df = self.get_df(mol=True)
		
		df.rename(inplace=True, columns={name_col:'_Name', 'mol':'ROMol'})

		logger.writing(out_path)
			
		from rdkit.Chem import PandasTools
		PandasTools.WriteSDF(df, out_path, "ROMol", "_Name", list(df.columns))

	def to_fragalysis(self, 
		out_path,

		*,
		method,
		ref_url,
		submitter_name,
		submitter_email,
		submitter_institution,
		metadata: bool = True,
	):

		"""Prepare an SDF for upload to the RHS of Fragalysis"""

		from .fragalysis import generate_header
		from rdkit.Chem import SDWriter, PandasTools

		name_col = '_Name'
		mol_col = 'ROMol'

		# get the dataframe of poses

		pose_df = self.get_df(mol=True, inspirations='fragalysis', reference='name', metadata=metadata)
		pose_df = pose_df.drop(columns=['id', 'path', 'compound', 'target', 'ref_pdb', 'original SMILES'], errors='ignore')

		pose_df.rename(inplace=True, columns=
			{'name':name_col, 
			 'mol':mol_col, 
			 'inspirations':'ref_mols',
			 'reference':'ref_pdb',
			 'smiles':'original SMILES',
			})

		# create the header molecule
		
		df_cols = set(pose_df.columns)

		header = generate_header(
			self[0],
			method=method,
			ref_url=ref_url,
			submitter_name=submitter_name,
			submitter_email=submitter_email,
			submitter_institution=submitter_institution,
			extras={'smiles':'smiles', 'ref_mols':'fragment inspirations'},
			metadata=metadata,
		)

		# return header

		header_cols = set(header.GetPropNames())

		# print(df_cols - header_cols)
		# print(header_cols - df_cols)
		# print(df_cols)

		# empty properties
		pose_df['generation_date'] = [None] * len(pose_df)
		pose_df['submitter_name'] = [None] * len(pose_df)
		pose_df['method'] = [None] * len(pose_df)
		pose_df['submitter_email'] = [None] * len(pose_df)
		pose_df['ref_url'] = [None] * len(pose_df)

		fields = []

		logger.writing(out_path)

		with open(out_path, 'w') as sdfh:
			with SDWriter(sdfh) as w:
				w.write(header)
			PandasTools.WriteSDF(pose_df, sdfh, mol_col, name_col, set(pose_df.columns))

		return pose_df

	def draw(self):
		"""Render this pose set with Py3Dmol"""
		
		from molparse.rdkit import draw_mols

		mols = [p.mol for p in self]

		return draw_mols(mols)

	def grid(self):
		"""Draw a grid of all contained molecules"""
		from molparse.rdkit import draw_grid

		data = [(p.name, p.compound.mol) for p in self]

		mols = [d[1] for d in data]
		labels = [d[0] for d in data]

		return draw_grid(mols, labels=labels)

	### OTHER

	def summary(self):
		"""Print a summary of this pose set"""
		logger.header('PoseSet()')
		logger.var('#poses', len(self))
		logger.var('#compounds', self.num_compounds)
		logger.var('tags', self.tags)

	### DUNDERS

	def __repr__(self) -> str:
		return f'{mcol.bold}{mcol.underline}subset(P x {len(self)}){mcol.unbold}{mcol.ununderline}'

	def __len__(self) -> int:
		return len(self.indices)

	def __iter__(self):
		return iter(self.db.get_pose(table=self.table, id=i) for i in self.indices)

	def __getitem__(self, key) -> Pose:
		try:
			index = self.indices[key]
		except IndexError:
			logger.exception(f'list index out of range: {key=} for {self}')
			raise
		return self.db.get_pose(table=self.table, id=index)
