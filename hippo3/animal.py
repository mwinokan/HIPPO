
import pandas as pd

from .cset import CompoundSet
from .pset import PoseSet
from .tags import TagSet
from .rset import ReactionSet

from .db import Database
from pathlib import Path

from tqdm import tqdm

from .tools import inchikey_from_smiles, sanitise_smiles

from mlog import setup_logger
logger = setup_logger('HIPPO')

class HIPPO:
		
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
		self._db = Database(self.db_path)

		self._compounds = CompoundSet(self.db, 'compound')
		self._poses = PoseSet(self.db, 'pose')
		self._tags = TagSet(self.db, 'tag')
		self._reactions = ReactionSet(self.db, 'reaction')

		logger.success(f"Initialised animal {self}")
		
	### FACTORIES

	### PROPERTIES

	@property
	def name(self):
		return self._name

	# @property
	# def target_name(self):
	# 	return self._target_name

	@property
	def db_path(self):
		return self._db_path

	@property
	def db(self):
		return self._db

	@property
	def compounds(self):
		return self._compounds

	@property
	def poses(self):
		return self._poses

	@property
	def reactions(self):
		return self._reactions

	@property
	def tags(self):
		return self._tags
	
	### MAIN METHODS

	def add_hits(self, target_name, metadata_csv, aligned_directory, skip=None, debug=False):

		""" Pass in a generator to the aligned directory from a Fragalysis download"""

		import molparse as mp
		from rdkit.Chem import PandasTools
		from .tools import remove_other_ligands

		# create the target
		target_id = self.db.insert_target(name=target_name)

		if not target_id:
			target_id = self.db.get_target_id(name=target_name)

		skip = skip or []

		logger.var('aligned_directory',aligned_directory)

		count_directories_tried = 0
		count_compound_registered = 0
		count_poses_registered = 0

		meta_df = pd.read_csv(metadata_csv)
		generated_tag_cols = ['ConformerSites', 'CanonSites', 'CrystalformSites', 'Quatassemblies', 'Crystalforms']
		curated_tag_cols = [c for c in meta_df.columns if c not in ['Code', 'Long code', 'Compound code', 'Smiles']+generated_tag_cols]

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
			df = PandasTools.LoadSDF(sdfs[0])

			# extract fields
			observation_shortname = path.name.replace('.sdf','')
			observation_longname = df.ID[0]
			mol = df.ROMol[0]
			lig_res_number = int(observation_longname.split('-')[1].split('_')[2])

			# create the single ligand bound pdb

			sys = mp.parse(pdbs[0], verbosity=debug)

			sys = remove_other_ligands(sys, lig_res_number)

			pose_path = str(pdbs[0].resolve()).replace('.pdb','_hippo.pdb')
			smiles = mp.rdkit.mol_to_smiles(mol)

			smiles = sanitise_smiles(smiles, verbosity=debug)

			mp.write(pose_path, sys, shift_name=True, verbosity=debug)

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

			metadata = {}
			for tag in generated_tag_cols:
				metadata[tag] = meta_row[tag].values[0]

			pose_id = self.db.insert_pose(
				compound=compound,
				name=observation_shortname,
				target=target_id,
				path=pose_path,
				tags=tags,
				metadata=metadata,
			)

			if pose_id:
				count_poses_registered += 1

		logger.var('#directories parsed', count_directories_tried)
		logger.var('#compounds registered', count_compound_registered)
		logger.var('#poses registered', count_poses_registered)

		return meta_df

	# def calculate_fingerprints(self, poses):

	# 	# group by reference ID??? Use an SQL select query to construct a PoseSubset with the groups

	# 	...



	### PLOTTING

	def plot_tag_statistics(self, *args, **kwargs):
		from .plotting import plot_tag_statistics
		return plot_tag_statistics(self, *args, **kwargs)

	def plot_compound_property(self, prop, **kwargs): 
		from .plotting import plot_compound_property
		return plot_compound_property(self, prop, **kwargs)

	def plot_pose_property(self, prop, **kwargs): 
		from .plotting import plot_pose_property
		return plot_pose_property(self, prop, **kwargs)

	def plot_interaction_punchcard(self, poses, subtitle, opacity=1.0, **kwargs):
		from .plotting import plot_interaction_punchcard
		return plot_interaction_punchcard(self, poses=poses, subtitle=subtitle, opacity=opacity, **kwargs)

	### DUNDERS

	def __repr__(self):
		return f'HIPPO("{self.name}")'
