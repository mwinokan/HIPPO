
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
		"""Returns the aliases of child poses"""
		return self.aliases

	@property
	def aliases(self):
		"""Returns the aliases of child poses"""
		result = self.db.select(table=self.table, query='pose_alias', multiple=True)
		return [q for q, in result]

	@property
	def inchikeys(self):
		"""Returns the inchikeys of child poses"""
		result = self.db.select(table=self.table, query='pose_inchikey', multiple=True)
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

	def draw(self, max_draw=100):
		if len(self) <= max_draw:
			self[:].draw()
		else:
			logger.warning(f"Too many poses: {len(self)} > {max_draw=}. Increase max_draw or use animal.poses[:].draw()")		

	def summary(self):
		"""Print a summary of this pose set"""
		logger.header('PoseTable()')
		logger.var('#poses', len(self))
		logger.var('tags', self.tags)

	### DUNDERS

	def __call__(self, tag=None):
		if tag:
			return self.get_by_tag(tag)

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
				pose = self.db.get_pose(alias=key)
				if not pose:
					pose = self.db.get_pose(inchikey=key)
				return pose

			case key if isinstance(key, list) or isinstance(key, tuple) or isinstance(key, set):

				indices = []
				for i in key:
					if isinstance(i,int):
						index = i
					elif isinstance(i,str):
						index = self.db.get_pose_id(alias=i)
						if not index:
							index = self.db.get_pose_id(inchikey=i)
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

		assert all(isinstance(i, int) for i in indices)

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
		"""Returns the aliases of poses in this set"""
		return self.aliases

	@property
	def aliases(self):
		"""Returns the aliases of child poses"""
		return [self.db.select_where(table=self.table, query='pose_alias', key='id', value=i, multiple=False)[0] for i in self.indices]

	@property
	def inchikeys(self):
		"""Returns the inchikeys of child poses"""
		return [self.db.select_where(table=self.table, query='pose_inchikey', key='id', value=i, multiple=False)[0] for i in self.indices]

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
	def mols(self):
		"""Get the rdkit Molecules contained in this set"""
		return [p.mol for p in self]

	@property
	def num_compounds(self):
		"""Count the compounds associated to this set of poses"""
		return len(self.compounds)

	@property
	def df(self):
		"""Get a DataFrame of the poses in this set"""
		return self.get_df(mol=True)

	@property
	def reference_ids(self):
		values = self.db.select_where(table='pose', query='DISTINCT pose_reference', key=f'pose_reference IS NOT NULL and pose_id in {tuple(self.ids)}', value=None, multiple=True)
		return set(v for v, in values)

	# @property
	# def inspirations(self):
	# 	return self._inspirations

	@property
	def inspiration_sets(self):
		sets = []
		for id in self.ids:
			insp_ids = self.db.select_where(query='inspiration_original', table='inspiration', key='derivative', value=id, multiple=True, sort='inspiration_original')
			insp_ids = set(v for v, in insp_ids)
			# print(insp_ids)
			if insp_ids not in sets:
				sets.append(insp_ids)
		return sets

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

	def get_by_reference(self, ref_id):
		"""Get poses with a certain reference id"""
		values = self.db.select_where(table='pose', query='pose_id', key='reference', value=ref_id, multiple=True)
		if not values:
			return None
		return PoseSet(self.db, [v for v, in values])

	def filter(self, function, inverse=False):
		"""Filter this poseset by selecting members where function(pose) is truthy"""
		
		ids = set()
		for pose in self:
			value = function(pose)
			# logger.debug(f'{pose=} {value=}')
			if value and not inverse:
				ids.add(pose.id)
			elif not value and inverse:
				ids.add(pose.id)

		return PoseSet(self.db, ids)

	### BULK SETTING

	@property
	def reference(self):
		raise NotImplementedError
	
	@reference.setter
	def reference(self, r):
		"""Bulk set the references for poses in this set"""
		if not isinstance(r, int):
			assert r._table == 'pose'
			r = r.id

		for i in self.indices:
			self.db.update(table='pose', id=i, key='pose_reference', value=r, commit=False)

		self.db.commit()
	
	def add_tag(self, t):
		"""Add this tag to every member of the set"""

		assert isinstance(t, str)

		for i in self.indices:
			self.db.insert_tag(name=t, pose=i, commit=False)	

		logger.info(f'Tagged {self} w/ "{t}"')		

		self.db.commit()

	def append_to_metadata(self, key, value):
		"""Create or append to a list-like value with given key for each pose in this set"""
		for id in self.indices:
			metadata = self.db.get_metadata(table='pose', id=id)
			metadata.append(key, value)

	### SPLITTING

	def split_by_reference(self):
		sets = {}
		for ref_id in self.reference_ids:
			sets[ref_id] = self.get_by_reference(ref_id)
		return sets

	def split_by_inspirations(self, single_set=False):
		
		sets = {}

		for pose in self:

			insp_ids = tuple(pose.get_inspiration_ids())

			if insp_ids not in sets:
				sets[insp_ids] = PoseSet(self.db, [pose.id])
			else:
				sets[insp_ids]._indices.append(pose.id)

		if single_set:
			logger.var('#unique inspiration combinations', len(sets))
			sets = PoseSet(self.db, sum([s.ids for s in sets.values()], []))

		return sets


	### EXPORTING

	def write_sdf(self, out_path, name_col='name'):
		"""Write an SDF"""

		from pathlib import Path

		df = self.get_df(mol=True)

		df.rename(inplace=True, columns={name_col:'_Name', 'mol':'ROMol'})

		logger.writing(out_path)
			
		from rdkit.Chem import PandasTools
		PandasTools.WriteSDF(df, out_path, "ROMol", "_Name", list(df.columns))

		# keep record of export
		value = str(Path(out_path).resolve())
		self.db.remove_metadata_list_item(table='pose', key='exports', value=value)
		self.append_to_metadata(key='exports', value=value)

	def to_fragalysis(self, 
		out_path,
		*,
		method,
		ref_url,
		submitter_name,
		submitter_email,
		submitter_institution,
		metadata: bool = True,
		sort_by: str | None = None,
		sort_reverse: bool = False,
		generate_pdbs: bool = False,
	):

		"""Prepare an SDF for upload to the RHS of Fragalysis"""

		from .fragalysis import generate_header
		from pathlib import Path
		from rdkit.Chem import SDWriter, PandasTools
		assert out_path.endswith('.sdf')

		name_col = '_Name'
		mol_col = 'ROMol'

		# get the dataframe of poses

		pose_df = self.get_df(mol=True, inspirations='fragalysis', duplicate_name='original ID', reference='name', metadata=metadata)
		pose_df = pose_df.drop(columns=['path', 'compound', 'target', 'ref_pdb', 'original SMILES'], errors='ignore')

		pose_df.rename(inplace=True, columns=
			{'name':name_col, 
			 'id':'HIPPO ID', 
			 'mol':mol_col, 
			 'inspirations':'ref_mols',
			 'reference':'ref_pdb',
			 'smiles':'original SMILES',
			})


		if generate_pdbs:
			
			from zipfile import ZipFile

			# output subdirectory
			out_key = Path(out_path).name.removesuffix('.sdf')
			pdb_dir = Path(out_key)
			pdb_dir.mkdir(exist_ok=True)

			# create the zip archive
			with ZipFile(f'{out_key}_pdbs.zip', 'w') as z:
				
				# loop over poses
				for (i,row),pose in zip(pose_df.iterrows(), self):

					# filenames
					pdb_name = f"{out_key}_{row._Name}.pdb"
					pdb_path = pdb_dir / pdb_name
					pose_df.iloc[i]['ref_pdb'] = pdb_name

					# generate the PL-complex
					sys = pose.complex_system
					
					# write the PDB
					logger.writing(pdb_path)
					sys.write(pdb_path, verbosity=0)
					z.write(pdb_path)
			
			logger.writing(f'{out_key}_pdbs.zip')

		# create the header molecule
		
		df_cols = set(pose_df.columns)

		header = generate_header(
			self[0],
			method=method,
			ref_url=ref_url,
			submitter_name=submitter_name,
			submitter_email=submitter_email,
			submitter_institution=submitter_institution,
			extras={'smiles':'smiles', 'ref_mols':'fragment inspirations', 'original ID':'original ID'},
			metadata=metadata,
		)

		header_cols = set(header.GetPropNames())

		# empty properties
		pose_df['generation_date'] = [None] * len(pose_df)
		pose_df['submitter_name'] = [None] * len(pose_df)
		pose_df['method'] = [None] * len(pose_df)
		pose_df['submitter_email'] = [None] * len(pose_df)
		pose_df['ref_url'] = [None] * len(pose_df)

		if sort_by:
			pose_df = pose_df.sort_values(by=sort_by, ascending=not sort_reverse)

		fields = []

		logger.writing(out_path)

		with open(out_path, 'w') as sdfh:
			with SDWriter(sdfh) as w:
				w.write(header)
			PandasTools.WriteSDF(pose_df, sdfh, mol_col, name_col, set(pose_df.columns))

		# keep record of export
		value = str(Path(out_path).resolve())
		self.db.remove_metadata_list_item(table='pose', key='exports', value=value)
		self.append_to_metadata(key='exports', value=value)

		return pose_df

	def to_pymol(self, prefix=None):
		"""Group the poses by reference protein and inspirations and output relevant PDBs and SDFs."""

		commands = []

		prefix = prefix or ''
		if prefix:
			prefix = f'{prefix}_'

		from pathlib import Path

		for i,(ref_id, poses) in enumerate(self.split_by_reference().items()):

			ref_pose = self.db.get_pose(id=ref_id)
			ref_name = ref_pose.name or ref_id
		
			# create the subdirectory
			ref_dir = Path(f'{prefix}ref_{ref_name}')
			logger.writing(ref_dir)
			ref_dir.mkdir(parents=True, exist_ok=True)
			
			# write the reference protein
			ref_pdb = ref_dir / f'ref_{ref_name}.pdb'
			ref_pose.protein_system.write(ref_pdb, verbosity=0)

			# color the reference:
			commands.append(f'load {ref_pdb.resolve()}')
			commands.append('hide')
			commands.append('show lines')
			commands.append('show surface')
			commands.append('util.cbaw')
			commands.append('set surface_color, white')
			commands.append('set transparency,  0.4')

			for j,(insp_ids, poses) in enumerate(poses.split_by_inspirations().items()):

				inspirations = PoseSet(self.db, insp_ids)
				insp_names = "-".join(inspirations.names)

				# create the subdirectory
				insp_dir = ref_dir / insp_names
				insp_dir.mkdir(parents=True, exist_ok=True)

				# write the inspirations
				insp_sdf = insp_dir / f'{insp_names}_frags.sdf'
				inspirations.write_sdf(insp_sdf)

				commands.append(f"load {insp_sdf.resolve()}")
				commands.append(f"set all_states, on, {insp_sdf.name.removesuffix('.sdf')}")
				commands.append(f"util.rainbow \"{insp_sdf.name.removesuffix('.sdf')}\"")

				# write the poses
				pose_sdf = insp_dir / f'{insp_names}_derivatives.sdf'
				poses.write_sdf(pose_sdf)
				
				commands.append(f"load {pose_sdf.resolve()}")
				commands.append(f'util.cbaw "{pose_sdf.name.removesuffix(".sdf")}"')

				if j > 0:
					commands.append(f"disable \"{insp_sdf.name.removesuffix('.sdf')}\"")
					commands.append(f'disable "{pose_sdf.name.removesuffix(".sdf")}"')

		return '; '.join(commands)

	### OUTPUT

	def interactive(self, method=None, print_name=True, **kwargs):
		"""Creates a ipywidget to interactively navigate this PoseSet."""

		from ipywidgets import interactive, BoundedIntText, Checkbox, interactive_output, HBox, GridBox, Layout, VBox
		from IPython.display import display
		from pprint import pprint

		if method:
			def widget(i):
				pose = self[i]
				if print_name:
					print(repr(pose))
				value = getattr(pose, method)(**kwargs)
				if value:
					display(value)

			return interactive(widget, i=
				BoundedIntText(
					value=0,
					min=0,
					max=len(self)-1,
					step=1,
					description='Pose:',
					disabled=False
				))

		else:

			a = BoundedIntText(
					value=0,
					min=0,
					max=len(self)-1,
					step=1,
					description=f'Pose (/{len(self)}):',
					disabled=False,
				)

			b = Checkbox(description='Name', value=True)
			c = Checkbox(description='Summary', value=False)
			d = Checkbox(description='2D', value=False)
			e = Checkbox(description='3D', value=True)
			f = Checkbox(description='Metadata', value=False)

			ui = GridBox([b, c, d, e, f], layout=Layout(grid_template_columns="repeat(5, 100px)"))
			ui = VBox([a, ui])
			
			def widget(i, name=True, summary=True, grid=True, draw=True, metadata=True):
				pose = self[i]
				if name:
					print(repr(pose))

				if summary: pose.summary(metadata=False)
				if grid: pose.grid()
				if draw: pose.draw()
				if metadata:
					logger.title('Metadata:')
					pprint(pose.metadata)

			out = interactive_output(widget, {'i': a, 'name': b, 'summary': c, 'grid':d, 'draw':e, 'metadata':f})

			display(ui, out)

	def summary(self):
		"""Print a summary of this pose set"""
		logger.header('PoseSet()')
		logger.var('#poses', len(self))
		logger.var('#compounds', self.num_compounds)
		logger.var('tags', self.tags)

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
