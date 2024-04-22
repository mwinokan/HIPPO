
import pandas as pd
# import numpy as np

from .cset import CompoundTable
from .pset import PoseTable, PoseSet
from .tags import TagTable
from .rset import ReactionTable
from .compound import Compound
from .reaction import Reaction
from .target import Target
from .pose import Pose

from .db import Database
from pathlib import Path

from tqdm import tqdm

from .tools import inchikey_from_smiles, sanitise_smiles, SanitisationError

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
		self._reactions = ReactionTable(self.db, 'reaction')

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
			
			# load the SDF
			df = PandasTools.LoadSDF(str(sdfs[0]), molColName='ROMol', idName='ID', strictParsing=True)

			# extract fields
			observation_shortname = path.name.replace('.sdf','')
			observation_longname = df.ID[0]
			mol = df.ROMol[0]
			lig_res_number = int(observation_longname.split('-')[1].split('_')[2])
			lig_chain = observation_longname.split('-')[1].split('_')[1]

			# smiles			
			smiles = mp.rdkit.mol_to_smiles(mol)
			smiles = sanitise_smiles(smiles, verbosity=debug)

			# parse PDB
			sys = mp.parse(pdbs[0], verbosity=debug)

			if len(sys.residues['LIG']) > 1 or any(r.contains_alternative_sites for r in sys.residues['LIG']):
				# create the single ligand bound pdb
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
				alias=observation_shortname,
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
		check_pose_RMSD: bool = False,
		pose_RMSD_tolerance: float = 1.0,
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

			if not comp:
				logger.error(f'Could not register compound {i=}')
				continue

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
				alias=name, 
				compound=comp, 
				target=target, 
				path=mol_path, 
				metadata=metadata, 
				inspirations=inspirations, 
				tags=tags, 
				reference=reference,
				check_RMSD=check_pose_RMSD,
				RMSD_tolerance=pose_RMSD_tolerance,
			)

			if stop_after and i+1 >= stop_after:
				break

		logger.success(f'Loaded compounds from {sdf_path}')

	def add_syndirella_elabs(self,
		df_path: str | Path,
		*,
		inspiration_map: dict | None = None,
		base_only: bool = False,
		tags: None | list = None,
		reaction_yield_map: dict | None = None,
	):

		tags = tags or []
		
		if isinstance(df_path, str):
			df_path = Path(df_path)
		df = pd.read_pickle(df_path)

		# work out number of reaction steps
		n_steps = max([int(s.split('_')[0]) for s in df.columns if '_product_smiles' in s])

		if f'{n_steps}_num_atom_diff' not in df.columns:
			logger.error(df_path)
			print(df.columns)
			raise NotImplementedError

		base_id = None

		if base_only:
			generator = df.iterrows()
		else:
			generator = tqdm(df.iterrows())

		n_comps = len(self.compounds)
		n_poses = len(self.poses)

		for i,row in generator:

			path_to_mol = row.path_to_mol

			# skip pose-less entries
			if not path_to_mol:
				continue
			
			path_to_mol = path_to_mol.replace('//','/')

			this_row_is_a_base = 'base' in row[f'{n_steps}_product_name']

			if base_only and not this_row_is_a_base:
				continue

			if base_only and base_id and not this_row_is_a_base:
				break

			try:
				
				# loop over each reaction step
				for j in range(n_steps):
	
					j += 1
	
					reactants = []
	
					reactant_previous_product = row[f'{j}_r_previous_product']
	
					# reactant 1
					if reactant_previous_product == 1:
						reactant1_id = product.id
						reactants.append(reactant1_id)
					elif smiles := row[f'{j}_r1_smiles']:
						
						if not isinstance(smiles, str):
							raise InvalidRowError(f'non-string {j}_r1_smiles')
							
						reactant1_id = self.register_compound(smiles=smiles, commit=False, return_compound=False, tags=tags)
						reactants.append(reactant1_id)
				
					# reactant 2
					if reactant_previous_product == 2:
						reactant2_id = product.id
						reactants.append(reactant2_id)
					elif smiles := row[f'{j}_r2_smiles']:
						
						if not isinstance(smiles, str):
							raise InvalidRowError(f'non-string {j}_r2_smiles')
							
						reactant2_id = self.register_compound(smiles=smiles, commit=False, return_compound=False, tags=tags)
						reactants.append(reactant2_id)
				
					# product
					if smiles := row[f'{j}_product_smiles']:
						if not isinstance(smiles, str):
							raise InvalidRowError(f'non-string {j}_product_smiles')
	
						if j != n_steps:
							this_tags = ['intermediate'] + tags
							base = None
	
						elif this_row_is_a_base:
							this_tags = ['base'] + tags
							base = None
	
						else:
							this_tags = ['elab'] + tags
							base = base_id
	
						product = self.register_compound(smiles=smiles, tags=this_tags, commit=False, return_compound=True, base=base)
						
						if not base_id and j == n_steps and this_row_is_a_base:
							base_id = product.id
	
					# register the reaction
					if reaction_yield_map:
						product_yield = reaction_yield_map[row[f'{j}_reaction']]
					else:
						product_yield = 1.0
						
					self.register_reaction(reactants=reactants, product=product, type=row[f'{j}_reaction'], commit=False, product_yield=product_yield)

			except InvalidRowError as e:
				logger.error(f'Skipping invalid row {i=}: {e}')
				continue

			# pose metadata
			metadata = {
				'ddG':row['∆∆G'],
				'RMSD':row['comRMSD'],
			}

			# inspirations
			inspirations = []
			for inspiration in row.regarded:
				if inspiration_map:
					inspiration = inspiration_map[inspiration]
				else:
					# this is really expensive
					inspiration = self.poses[inspiration]
				if inspiration:
					inspirations.append(inspiration)
			
			# register the pose
			self.register_pose(
				compound=product, 
				target='A71EV2A', 
				path=path_to_mol, 
				inspirations=inspirations, 
				metadata=metadata, 
				tags=this_tags, 
				commit=False, 
				return_pose=False,
				overwrite_metadata=True,
				warn_duplicate=False,
			)

		self.db.commit()

		n_comps = len(self.compounds) - n_comps
		n_poses = len(self.poses) - n_poses

		if n_comps:
			if not base_only: 
				logger.success(f'Loaded {n_comps} new compounds from {df_path.name}')
		else:
			logger.warning(f'Loaded {n_comps} new compounds from {df_path.name}')
		if n_poses:
			if not base_only: 
				logger.success(f'Loaded {n_poses} new poses from {df_path.name}')
		else:
			logger.warning(f'Loaded {n_poses} new poses from {df_path.name}')

		return base_id

	### SINGLE INSERTION

	def register_compound(self, 
		*,
		smiles: str, 
		base: Compound | int | None = None, 
		tags: None | list = None, 
		metadata: None | dict = None,
		return_compound: bool = True,
		commit: bool = True,
		alias: str | None = None,
	) -> Compound:

		"""Use a smiles string to add a compound to the database. If it already exists return the compound"""

		assert smiles
		assert isinstance(smiles, str), f'Non-string {smiles=}'
		
		try:
			smiles = sanitise_smiles(smiles, sanitisation_failed='error')
		except SanitisationError:
			logger.error(f'Could not sanitise {smiles=}')
			return None
		except AssertionError:
			logger.error(f'Could not sanitise {smiles=}')
			return None

		inchikey = inchikey_from_smiles(smiles)

		compound_id = self.db.insert_compound(
			smiles=smiles, base=base, inchikey=inchikey, tags=tags, 
			metadata=metadata, warn_duplicate=False, commit=commit, alias=alias)

		if return_compound or metadata or alias or tags:
			if not compound_id:
				compound = self.compounds[inchikey]
			else:
				compound = self.compounds[compound_id]

			if metadata:
				compound.metadata.update(metadata)

			if alias:
				compound.alias = alias

			if tags:
				for tag in tags:
					compound.tags.add(tag)

			if return_compound:
				return compound
			else:
				return compound.id

		else:
			if not compound_id:
				compound_id = self.db.get_compound_id(inchikey=inchikey)
			return compound_id

	def register_reaction(self,
		*,
		type: str,
		product: Compound,
		reactants: list[Compound | int],
		commit: bool = True,
		product_yield: float = 1.0,
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

		reaction_id = self.db.insert_reaction(type=type, product=product, commit=commit, product_yield=product_yield)

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
		inchikey: str | None = None,
		alias: str | None = None,
		reference: int | None = None,
		tags: None | list = None,
		metadata: None | dict = None,
		inspirations: None | list = None,
		return_pose: bool = True,
		energy_score: float | None = None,
		distance_score: float | None = None,
		commit: bool = True,
		overwrite_metadata: bool = True,
		warn_duplicate: bool = True,
		check_RMSD: bool = False,
		RMSD_tolerance: float = 1.0,
	) -> Pose:

		if isinstance(compound, int):
			compound_id = compound
		else:
			compound_id = compound.id

		if check_RMSD:
			
			# check if the compound has existing poses
			other_pose_ids = self.db.select_id_where(table='pose', key='compound', value=compound_id, none='quiet', multiple=True)
			
			if other_pose_ids:
				other_poses = PoseSet(self.db, [i for i, in other_pose_ids])

				from molparse.rdkit import draw_mols, draw_flat
				from rdkit.Chem import MolFromMolFile
				from numpy.linalg import norm
				from numpy import array
				mol = MolFromMolFile(str(path.resolve()))
				
				c1 = mol.GetConformer()
				atoms1 = [a for a in mol.GetAtoms()]
				symbols1 = [a.GetSymbol() for a in atoms1]
				positions1 = [c1.GetAtomPosition(i) for i,_ in enumerate(atoms1)]

				for pose in other_poses:
					c2 = pose.mol.GetConformer()
					atoms2 = [a for a in pose.mol.GetAtoms()]
					symbols2 = [a.GetSymbol() for a in atoms2]
					positions2 = [c2.GetAtomPosition(i) for i,_ in enumerate(atoms2)]
					
					for s1,p1 in zip(symbols1, positions1):
						for s2,p2 in zip(symbols2, positions2):
							if s2 != s1:
								continue
							if norm(array(p2-p1)) <= RMSD_tolerance:
								# this atom (1) is within tolerance
								break
						else:
							# this atom (1) is outside of tolerance
							break
					else:
						# all atoms within tolerance --> too similar
						logger.warning(f'Found similar {pose=}')
						if return_pose:
							return pose
						else:
							return pose.id
					
		pose_id = self.db.insert_pose(
			compound=compound, 
			inchikey=inchikey, 
			alias=alias, 
			target=target, 
			path=path, 
			tags=tags, 
			metadata=metadata, 
			reference=reference, 
			warn_duplicate=warn_duplicate, 
			commit=commit, 
			energy_score=energy_score, 
			distance_score=distance_score,
		)
		
		if not pose_id:
			if isinstance(path, Path):
				path = path.resolve()
			pose_id, = self.db.select_where(table='pose', query='pose_id', key='path', value=str(path))

		if not pose_id:
			logger.var('compound', compound)
			logger.var('inchikey', inchikey)
			logger.var('alias', alias)
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

class InvalidRowError(Exception):
	...
