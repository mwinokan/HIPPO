
import pandas as pd

from .cset import CompoundTable
from .pset import PoseTable
from .tags import TagTable
from .rset import ReactionSet
from .compound import Compound
from .reaction import Reaction
from .target import Target
from .pose import Pose

from .db import Database
from pathlib import Path

from tqdm import tqdm

from .tools import inchikey_from_smiles, sanitise_smiles

from mlog import setup_logger
logger = setup_logger('HIPPO')

class HIPPO:

	"""The HIPPO 'animal' class. Instantiating a HIPPO object will create or link a HIPPO database.

	::
		
		from hippo import HIPPO
		animal = HIPPO(project_name, db_path)

	* *project_name* give this hippo a name
	* *db_path* path where the database will be stored


	"""
		
	def __init__(self, 
		name: str,  
		db_path: str | Path,
	):

		logger.header('Creating HIPPO animal')

		self._name = name
		# self._target_name = target

		logger.var('name', name, dict(color='arg'))
		# logger.var('target', target, dict(color='arg'))

		if not isinstance(db_path, Path):
			db_path = Path(db_path)
		
		logger.var('db_path', db_path, dict(color='file'))

		self._db_path = db_path
		self._db = Database(self.db_path, animal=self)

		self._compounds = CompoundTable(self.db, 'compound')
		self._poses = PoseTable(self.db, 'pose')
		self._tags = TagTable(self.db, 'tag')
		self._reactions = ReactionSet(self.db, 'reaction')

		logger.success(f"Initialised animal {self}")
		
	### FACTORIES

	### PROPERTIES

	@property
	def name(self) -> str:
		"""Returns the project name"""
		return self._name

	@property
	def db_path(self) -> str:
		"""Returns the database path"""
		return self._db_path

	@property
	def db(self):
		"""Returns the Database object"""
		return self._db

	@property
	def compounds(self):
		"""Returns a :doc:`CompoundTable <compounds>` object, which interfaces to the compound table in the database"""
		return self._compounds

	@property
	def poses(self):
		"""Returns a :doc:`PoseTable <poses>` object, which interfaces to the pose table in the database"""
		return self._poses

	@property
	def reactions(self):
		"""Returns a :doc:`ReactionSet <reactions>` object, which interfaces to the reaction table in the database"""
		return self._reactions

	@property
	def tags(self):
		"""Returns a :doc:`TagTable <metadata>` object, which interfaces to the tag table in the database"""
		return self._tags
	
	@property
	def num_compounds(self):
		"""Returns the number of compounds"""
		return len(self.compounds)

	@property
	def num_poses(self):
		"""Returns the number of poses"""
		return len(self.poses)

	@property
	def num_reactions(self):
		"""Returns the number of reactions"""
		return len(self.reactions)

	@property
	def num_tags(self):
		"""Returns the number of tags"""
		return len(self.tags.unique)

	### BULK INSERTION

	def add_hits(self, 
		target_name: str, 
		metadata_csv: str | Path, 
		aligned_directory: str | Path, 
		skip: list | None = None, 
		debug: bool = False, 
	) -> pd.DataFrame:

		"""Load in crystallographic hits downloaded from Fragalysis.

		* *target_name* Name of this protein target
		* *metadata_csv* Path to the metadata.csv from the Fragalysis download
		* *aligned_directory* Path to the aligned_files directory from the Fragalysis download
		* *skip* optional list of observation names to skip

		Returns: a DataFrame of metadata

		"""

		import molparse as mp
		from rdkit.Chem import PandasTools
		from .tools import remove_other_ligands

		if not isinstance(aligned_directory, Path):
			aligned_directory = Path(aligned_directory)

		# create the target
		target = self.register_target(name=target_name)

		skip = skip or []

		logger.var('aligned_directory',aligned_directory)

		count_directories_tried = 0
		count_compound_registered = 0
		count_poses_registered = 0

		meta_df = pd.read_csv(metadata_csv)
		generated_tag_cols = ['ConformerSites', 'CanonSites', 'CrystalformSites', 'Quatassemblies', 'Crystalforms']
		curated_tag_cols = [c for c in meta_df.columns if c not in ['Code', 'Long code', 'Compound code', 'Smiles', 'Downloaded']+generated_tag_cols]

		for path in tqdm(aligned_directory.iterdir()):

			if not path.is_dir():
				continue

			if path.name in skip:
				continue

			count_directories_tried += 1

			sdfs = list(path.glob('*x?????.sdf'))
			pdbs = list(path.glob('*x?????.pdb'))

			assert len(sdfs) == 1, (path, sdfs)
			assert len(pdbs) == 1, (path, pdbs)

			# logger.var('sdfs[0]',sdfs[0])
			# logger.var('PandasTools', PandasTools.__file__)
			
			# load the SDF
			df = PandasTools.LoadSDF(str(sdfs[0]), molColName='ROMol', idName='ID', strictParsing=True)

			# extract fields
			observation_shortname = path.name.replace('.sdf','')
			observation_longname = df.ID[0]
			mol = df.ROMol[0]
			lig_res_number = int(observation_longname.split('-')[1].split('_')[2])
			lig_chain = observation_longname.split('-')[1].split('_')[1]
			
			# logger.var('observation_longname',observation_longname)

			# smiles			
			smiles = mp.rdkit.mol_to_smiles(mol)
			smiles = sanitise_smiles(smiles, verbosity=debug)

			# parse PDB
			sys = mp.parse(pdbs[0], verbosity=debug)

			if len(sys.residues['LIG']) > 1 or any(r.contains_alternative_sites for r in sys.residues['LIG']):
				# create the single ligand bound pdb
				# logger.debug(str(pdbs[0].resolve()))
				sys = remove_other_ligands(sys, lig_res_number, lig_chain)
				sys.prune_alternative_sites('A', verbosity=0)
				pose_path = str(pdbs[0].resolve()).replace('.pdb', '_hippo.pdb')
				mp.write(pose_path, sys, shift_name=True, verbosity=debug)
			else:
				pose_path = str(pdbs[0].resolve())
			
			# create the molecule / pose
			compound_id = self.db.insert_compound(
				# name=crystal_name, 
				smiles=smiles, 
				tags=['hits'],
				warn_duplicate=debug,
			)

			if not compound_id:

				inchikey = inchikey_from_smiles(smiles)
				compound = self.compounds[inchikey]

				if not compound:
					logger.error('Compound exists in database but could not be found by inchikey')
					logger.var('smiles',smiles)
					logger.var('inchikey',inchikey)
					logger.var('observation_shortname',observation_shortname)
					raise Exception
			
			else:
				count_compound_registered += 1
				compound = self.compounds[compound_id]

			# metadata

			meta_row = meta_df[meta_df['Code']==observation_shortname]

			tags = ['hits']
			for tag in curated_tag_cols:
				if meta_row[tag].values[0]:
					tags.append(tag)

			metadata = {'observation_longname':observation_longname}
			for tag in generated_tag_cols:
				metadata[tag] = meta_row[tag].values[0]

			pose_id = self.db.insert_pose(
				compound=compound,
				name=observation_shortname,
				target=target.id,
				path=pose_path,
				tags=tags,
				metadata=metadata,
			)

			if pose_id:
				count_poses_registered += 1

		logger.var('#directories parsed', count_directories_tried)
		logger.var('#compounds registered', count_compound_registered)
		logger.var('#poses registered', count_poses_registered)

	def add_compounds(self, 
		target: str,
		sdf_path: str | Path, 
		*,
		reference: int | Pose | None = None,
		tags: None | list = None,
		output_directory: str | Path | None = None,
		mol_col='ROMol',
		name_col='ID',
		inspiration_col = 'ref_mols',
		inspiration_map: None | dict = None,
		skip_first=False,
		convert_floats=True,
		stop_after: int | None = None,
		skip_equal_dict: dict | None = None,
		skip_not_equal_dict: dict | None = None,
	):

		"""Add virtual hits from an SDF into the database.

		All columns except those specified via arguments mol_col and name_col are added to the Pose metadata.
		"""

		if not isinstance(sdf_path, Path):
			sdf_path = Path(sdf_path)

		logger.debug(f'{sdf_path=}')

		tags = tags or []

		from rdkit.Chem import PandasTools, MolToMolFile, MolFromMolFile
		from molparse.rdkit import mol_to_smiles, mol_to_pdb_block
		from numpy import isnan
		from pandas import read_pickle

		if sdf_path.name.endswith('.sdf'):
			df = PandasTools.LoadSDF(sdf_path)
		else:
			df = read_pickle(sdf_path)

		df_columns = list(df.columns)

		target = self.register_target(target).name

		assert mol_col in df_columns, f'{mol_col=} not in {df_columns}'
		assert name_col in df_columns, f'{name_col=} not in {df_columns}'
		assert inspiration_col in df_columns, f'{inspiration_col=} not in {df_columns}'

		df_columns.pop(df_columns.index(mol_col))
		df_columns.pop(df_columns.index(name_col))
		df_columns.pop(df_columns.index(inspiration_col))

		if not output_directory:
			import os
			output_directory = str(sdf_path.name).removesuffix('.sdf')
			logger.writing(f'Creating output directory {output_directory}')
			os.system(f'mkdir -p {output_directory}')

		output_directory = Path(output_directory)

		for i,row in tqdm(df.iterrows()):

			name = row[name_col].strip()

			if name == 'ver_1.2':
				logger.warning('Skipping Fragalysis header molecule')
				continue

			if skip_equal_dict and any(row[k] == v for k,v in skip_equal_dict.items()):
				continue

			if skip_not_equal_dict and any(row[k] != v for k,v in skip_not_equal_dict.items()):
				continue
			
			mol = row[mol_col]

			if not name:
				name = f'pose_{i}'

			mol_path = output_directory / f'{name}.mol'

			MolToMolFile(mol, mol_path)

			smiles = mol_to_smiles(mol)

			comp = self.register_compound(smiles=smiles, tags=tags)

			# inspirations

			inspirations = []

			insp_str = row[inspiration_col]

			if isinstance(insp_str, str):
				insp_str = insp_str.removeprefix('[')
				insp_str = insp_str.removesuffix(']')
				insp_str = insp_str.replace("'", "")
				generator = insp_str.split(',')

			else:
				generator = insp_str

			for insp in generator:
				insp = insp.strip()

				if inspiration_map:
					if isinstance(inspiration_map, dict) and insp in inspiration_map:
						pose = inspiration_map[insp]
						if pose:
							inspirations.append(pose)
					elif hasattr(inspiration_map, '__call__'):
						pose = inspiration_map(insp)
						if pose:
							inspirations.append(pose)
				else:
					pose = self.poses[insp]
					if pose:
						inspirations.append(pose)
					else:
						logger.error(f'Could not find inspiration pose {insp}')
						continue

			# metadata
			metadata = {}

			for col in df_columns:
				value = row[col]

				if not isinstance(col, str):
					if i == 0:
						logger.warning(f'Skipping metadata from column={col}.')
					continue

				if isinstance(value, float) and isnan(value):
					continue

				if convert_floats:
					try:
						value = float(value)
					except TypeError:
					# 	logger.warning(f'Could not convert column={col}.')
						pass
					except ValueError:
						pass

				if not (isinstance(value, str) or isinstance(value, float)):
					if i == 0:
						logger.warning(f'Skipping metadata from column={col}.')
					continue

				metadata[col] = value

			# print(metadata)

			pose = self.register_pose(
				name=name, 
				compound=comp, 
				target=target, 
				path=mol_path, 
				metadata=metadata, 
				inspirations=inspirations, 
				tags=tags, 
				reference=reference
			)

			if stop_after and i+1 >= stop_after:
				break

		logger.success(f'Loaded compounds from {sdf_path}')

	def add_syndirella_elabs(self,
		df_path: str | Path,
		*,
		inspiration_map: dict | None = None,

	):

		logger.warning('Assuming this is a three-step reaction')

		df = pd.read_pickle(df_path)

		for i,row in tqdm(df.iterrows()):

			# skip pose-less entries

			if not row.path_to_mol:
				continue

			if i > 100:
				break

			# loop over each reaction step
			for j in range(3):

				j += 1

				reactants = []

				reactant_previous_product = row[f'{j}_r_previous_product']

				# reactant 1
				if reactant_previous_product == 1:
					reactant1_id = product.id
					reactants.append(reactant1_id)
				elif smiles := row[f'{j}_r1_smiles']:
					reactant1_id = self.register_compound(smiles=smiles, commit=False, return_compound=False)
					reactants.append(reactant1_id)
			
				# reactant 2
				if reactant_previous_product == 2:
					reactant2_id = product.id
					reactants.append(reactant2_id)
				elif smiles := row[f'{j}_r2_smiles']:
					reactant2_id = self.register_compound(smiles=smiles, commit=False, return_compound=False)
					reactants.append(reactant2_id)
			
				# product
				if smiles := row[f'{j}_product_smiles']:
					if j != 3:
						tags = ['intermediate']
					elif row['3_num_atom_diff'] == 0:
						tags = ['base']
					else:
						tags = ['elab']
					product = self.register_compound(smiles=smiles, tags=tags, commit=False, return_compound=True)
			
				# register the reaction
				self.register_reaction(reactants=reactants, product=product, type=row[f'{j}_reaction'], commit=False)

			# pose metadata
			metadata = {
				'ddG':row['∆∆G'],
				'RMSD':row['comRMSD'],
			}

			# inspirations
			inspirations = []
			
			# register the pose
			self.register_pose(
				compound=product, 
				target='A71EV2A', 
				path=row.path_to_mol, 
				inspirations=inspirations, 
				metadata=metadata, 
				tags=tags, 
				commit=False, 
				return_pose=False,
				overwrite_metadata=True,
			)

		self.db.commit()
		logger.success('Committed changes')

	### SINGLE INSERTION

	def register_compound(self, 
		*,
		smiles: str, 
		base: Compound | int | None = None, 
		tags: None | list = None, 
		metadata: None | dict = None,
		return_compound: bool = True,
		commit: bool = True,
	) -> Compound:

		"""Use a smiles string to add a compound to the database. If it already exists return the compound"""

		smiles = sanitise_smiles(smiles)
		inchikey = inchikey_from_smiles(smiles)

		compound_id = self.db.insert_compound(
			smiles=smiles, base=base, inchikey=inchikey, tags=tags, 
			metadata=metadata, warn_duplicate=False, commit=commit)

		if return_compound or metadata:
			if not compound_id:
				compound = self.compounds[inchikey]
			else:
				compound = self.compounds[compound_id]

			if metadata:
				compound.metadata.update(metadata)

			if return_compound:
				return compound
			else:
				return compound.id

		else:
			if not compound_id:
				compound_id = self.db.get_compound_id(name=inchikey)
			return compound_id

	def register_reaction(self,
		*,
		type: str,
		product: Compound,
		reactants: list[Compound | int],
		commit: bool = True,
	) -> Reaction:

		### CHECK FOR DUPLICATES
		
		reactant_ids = set(v if isinstance(v, int) else v.id for v in reactants)
		
		for reaction in product.reactions:

			if reaction.type != type:
				continue

			if reaction.product != product:
				continue

			if reaction.reactant_ids != reactant_ids:
				continue

			return reaction

		reaction_id = self.db.insert_reaction(type=type, product=product, commit=commit)

		for reactant in reactants:
			self.db.insert_reactant(compound=reactant, reaction=reaction_id, commit=commit)

		return self.reactions[reaction_id]

	def register_target(self,
		name: str,
	) -> Target:
		
		target_id = self.db.insert_target(name=name)

		if not target_id:
			target_id = self.db.get_target_id(name=name)

		return self.db.get_target(id=target_id)

	def register_pose(self,
		*,
		compound: Compound | int,
		target: str,
		path: str,
		name: str | None = None,
		reference: int | None = None,
		tags: None | list = None,
		metadata: None | dict = None,
		inspirations: None | list = None,
		return_pose: bool = True,
		commit: bool = True,
		overwrite_metadata: bool = True,
	) -> Pose:

		if isinstance(compound, int):
			compound_id = compound
		else:
			compound_id = compound.id
		
		pose_id = self.db.insert_pose(
			compound=compound, name=name, target=target, path=path, 
			tags=tags, metadata=metadata, reference=reference, warn_duplicate=True, commit=commit)
		
		if not pose_id:
			if isinstance(path, Path):
				path = path.resolve()
			pose_id, = self.db.select_where(table='pose', query='pose_id', key='path', value=str(path))

		if not pose_id:
			logger.var('compound', compound)
			logger.var('name', name)
			logger.var('target', target)
			logger.var('path', path)
			logger.var('reference', reference)
			logger.var('tags', tags)
			logger.debug(f'{metadata=}')
			logger.debug(f'{inspirations=}')

			raise Exception

		if return_pose or (metadata and not overwrite_metadata):
			pose = self.poses[pose_id]

			if metadata:
				pose.metadata.update(metadata)
		
		else:
			pose = pose_id

		if overwrite_metadata:
			self.db.insert_metadata(table='pose', id=compound_id, payload=metadata, commit=commit)

		inspirations = inspirations or []
		for inspiration in inspirations:
			self.db.insert_inspiration(original=inspiration, derivative=pose, warn_duplicate=False, commit=commit)

		return pose

	### PLOTTING

	def plot_tag_statistics(self, *args, **kwargs):
		"""Plot an overview of the number of compounds and poses for each tag"""
		from .plotting import plot_tag_statistics
		return plot_tag_statistics(self, *args, **kwargs)

	def plot_compound_property(self, prop, **kwargs): 
		"""Plot an arbitrary compound property across the whole dataset"""
		from .plotting import plot_compound_property
		return plot_compound_property(self, prop, **kwargs)

	def plot_pose_property(self, prop, **kwargs): 
		"""Plot an arbitrary pose property across the whole dataset"""
		from .plotting import plot_pose_property
		return plot_pose_property(self, prop, **kwargs)

	def plot_interaction_punchcard(self, poses, subtitle, opacity=1.0, **kwargs):
		"""Plot an interaction punchcard for a set of poses"""
		from .plotting import plot_interaction_punchcard
		return plot_interaction_punchcard(self, poses=poses, subtitle=subtitle, opacity=opacity, **kwargs)

	def plot_residue_interactions(self, poses, residue_number, **kwargs):
		"""Plot an interaction punchcard for a set of poses"""
		from .plotting import plot_residue_interactions
		return plot_residue_interactions(self, poses=poses, residue_number=residue_number, **kwargs)

	### OTHER

	def summary(self):
		"""Print a summary of this HIPPO"""
		logger.header(self)
		logger.var('db_path', self.db_path)
		logger.var('#compounds', self.num_compounds)
		logger.var('#poses', self.num_poses)
		logger.var('#reactions', self.num_reactions)
		logger.var('#tags', self.num_tags)
		logger.var('tags', self.tags.unique)

	### DUNDERS

	def __repr__(self) -> str:
		"""Returns a command line representation"""
		return f'HIPPO("{self.name}")'
