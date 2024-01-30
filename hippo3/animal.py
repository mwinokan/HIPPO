
import pandas as pd

from .cset import CompoundSet

from .db import Database
from pathlib import Path

from mlog import setup_logger
logger = setup_logger('HIPPO', debug=True)

class HIPPO:
		
	def __init__(self, name, target, db_path=None):

		logger.header('Creating HIPPO animal')

		self._name = name
		self._target_name = target

		logger.var('name', name, dict(color='arg'))
		logger.var('target', target, dict(color='arg'))

		db_path = db_path or f'{name}.db'
		db_path = Path(db_path)
		
		logger.var('db_path', db_path, dict(color='file'))

		self._db_path = db_path
		self._db = Database(self.db_path)

		self._compounds = CompoundSet(self.db, 'compound')

		logger.success(f"Initialised animal {self}")
		
	### FACTORIES

	### PROPERTIES

	@property
	def name(self):
		return self._name

	@property
	def target_name(self):
		return self._target_name

	@property
	def db_path(self):
		return self._db_path

	@property
	def db(self):
		return self._db

	@property
	def compounds(self):
		return self._compounds
	
	### PUBLIC METHODS

	def add_hits(self, metadata_path, pdb_path, pdb_pattern='**/*-x????_??_bound.pdb', tags=None, overwrite=False):
				
		### checks
		
		if not overwrite and 'hits' in [s.name for s in self.compound_sets]:
			logger.error(f'CompoundSet "hits" already exists')
			return

		if not isinstance(metadata_path, Path):
			metadata_path = Path(metadata_path)

		if not isinstance(pdb_path, Path):
			pdb_path = Path(pdb_path)

		tags = tags or ['hits']

		### metadata

		logger.var('metadata_path',str(metadata_path),dict(color='file'))

		try:
			metadata_df = pd.read_csv(metadata_path)
		except FileNotFoundError as e:
			logger.exception(e)
			raise

		# logger.header('metadata columns:')
		# pprint(list(metadata_df.columns))

		logger.success(f'Parsed metadata CSV')

		### pdbs

		logger.var('pdb_path',str(pdb_path),dict(color='file'))
		logger.var('pdb_pattern',str(pdb_pattern),dict(color='file'))

		pdbs = list(pdb_path.glob(pdb_pattern))

		pdbs = sorted(pdbs)

		if len(pdbs) < 1:
			logger.error(f'Did not find any PDBs',fatal=True,code='HIPPO.add_hits.0')
		
		# from .io import compounds_from_bound_pdbs

		comp_set = CompoundSet.from_bound_pdbs(name='hits', 
			pdbs=pdbs, 
			metadata_df=metadata_df,
			pdb_pattern=pdb_pattern,
			tags=tags,
			animal=self,
		)

		# if 'hits' in self.compound_sets:
		# 	self._compound_sets['hits'] = comp_set
		# else:
		# 	self._compound_sets.append(comp_set)
		
		# self._hit_compounds = self._compound_sets['hits']

		# logger.success(f'Loaded {comp_set.num_compounds} compounds as "{comp_set.name}" ({self.hits.num_poses} poses)')

	### DUNDERS

	def __repr__(self):
		return f'HIPPO("{self.name}")'
