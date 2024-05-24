
import mout
import mcol
import sqlite3
from sqlite3 import Error
from pprint import pprint
import json

from .compound import Compound
from .pose import Pose
from .reaction import Reaction
from .quote import Quote
from .metadata import MetaData
from .target import Target
from .feature import Feature

from .tools import inchikey_from_smiles

from pathlib import Path

import logging
logger = logging.getLogger('HIPPO')

class Database:

	"""Wrapper to connect to the HIPPO sqlite database."""

	def __init__(self, path: Path, animal) -> None:

		assert isinstance(path, Path)
		
		logger.debug('hippo3.Database.__init__()')

		self._path = path
		self._connection = None
		self._cursor = None
		self._animal = animal

		logger.debug(f'Database.path = {self.path}')

		try:
			path = path.resolve(strict=True)

		except FileNotFoundError:
			# create a blank database
			self.connect()
			self.create_blank_db()

		else:
			# existing database
			self.connect()
		
	### PROPERTIES

	@property
	def path(self) -> Path:
		"""Returns the path to the database file"""
		return self._path

	@property
	def connection(self):
		"""Returns a sqlite3.connection"""
		if not self._connection:
			self.connect()
		return self._connection
	
	@property
	def cursor(self):
		"""Returns a sqlite3.cursor"""
		if not self._cursor:
			self._cursor = self.connection.cursor()
		return self._cursor

	@property
	def total_changes(self):
		"""Return the total number of database rows that have been modified, inserted, or deleted since the database connection was opened."""
		return self.connection.total_changes

	### PUBLIC METHODS / API CALLS

	def close(self) -> None:
		"""Close the connection"""
		logger.debug('hippo3.Database.close()')
		if self.connection:
			self.connection.close()
		mout.success(f'Closed connection to {mcol.file}{self.path}')

	### GENERAL SQL

	def connect(self):
		"""Connect to the database"""
		logger.debug('hippo3.Database.connect()')

		conn = None

		try:
			conn = sqlite3.connect(self.path)

			logger.debug(f'{sqlite3.version=}')

			conn.enable_load_extension(True)
			conn.load_extension('chemicalite')
			conn.enable_load_extension(False)

			logger.success(f'Database connected @ {mcol.file}{self.path}')
		
		except sqlite3.OperationalError as e:

			if 'cannot open shared object file' in str(e):
				logger.error('chemicalite package not installed correctly')
			else:
				logger.exception(e)
			raise

		except Error as e:
			logger.exception(e)
			raise
		
		self._connection = conn
		self._cursor = conn.cursor()

	def execute(self, sql, payload=None):
		"""Execute arbitrary SQL"""
		try:
			if payload:
				return self.cursor.execute(sql, payload)
			else:
				return self.cursor.execute(sql)
		except Error as e:
			# logger.exception(e)
			raise

	def commit(self):
		"""Commit the changes"""
		# logger.debug('db.commit...')
		self.connection.commit()
		# raise Exception

	### CREATE TABLES

	def create_blank_db(self):
		"""Create a blank database"""

		logger.out('Creating blank database...')
		self.create_table_compound()
		self.create_table_inspiration()
		self.create_table_reaction()
		self.create_table_reactant()
		self.create_table_pose()
		self.create_table_tag()
		self.create_table_quote()
		self.create_table_target()
		self.create_table_pattern_bfp()
		self.create_table_feature()
		# self.create_table_morgan_bfp()

	def create_table_compound(self):
		"""Create the compound table"""
		logger.debug('HIPPO.Database.create_table_compound()')

		sql = """CREATE TABLE compound(
			compound_id INTEGER PRIMARY KEY,
			compound_inchikey TEXT,
			compound_alias TEXT,
			compound_smiles TEXT,
			compound_base INTEGER,
			compound_mol MOL,
			compound_pattern_bfp bits(2048),
			compound_morgan_bfp bits(2048),
			compound_metadata TEXT,
			FOREIGN KEY (compound_base) REFERENCES compound(compound_id),
			CONSTRAINT UC_compound_inchikey UNIQUE (compound_inchikey)
			CONSTRAINT UC_compound_alias UNIQUE (compound_alias)
			CONSTRAINT UC_compound_smiles UNIQUE (compound_smiles)
		);
		"""
			# CONSTRAINT UC_compound_morgan_bfp UNIQUE (compound_morgan_bfp)
			# CONSTRAINT UC_compound_pattern_bfp UNIQUE (compound_pattern_bfp)

		self.execute(sql)

	def create_table_inspiration(self):
		"""Create the inspiration table"""
		logger.debug('HIPPO.Database.create_table_inspiration()')

		sql = """CREATE TABLE inspiration(
			inspiration_original INTEGER,
			inspiration_derivative INTEGER,
			FOREIGN KEY (inspiration_original) REFERENCES pose(pose_id),
			FOREIGN KEY (inspiration_derivative) REFERENCES pose(pose_id),
			CONSTRAINT UC_inspiration UNIQUE (inspiration_original, inspiration_derivative)
		);
		"""

		self.execute(sql)

	def create_table_reaction(self):
		"""Create the reaction table"""
		logger.debug('HIPPO.Database.create_table_reaction()')

		sql = """CREATE TABLE reaction(
			reaction_id INTEGER PRIMARY KEY,
			reaction_type TEXT,
			reaction_product INTEGER,
			reaction_product_yield REAL,
			FOREIGN KEY (reaction_product) REFERENCES compound(compound_id)
		);
		"""

		self.execute(sql)

	def create_table_reactant(self):
		"""Create the reactant table"""
		logger.debug('HIPPO.Database.create_table_reactant()')

		sql = """CREATE TABLE reactant(
			reactant_amount REAL,
			reactant_reaction INTEGER,
			reactant_compound INTEGER,
			FOREIGN KEY (reactant_reaction) REFERENCES reaction(reaction_id),
			FOREIGN KEY (reactant_compound) REFERENCES compound(compound_id),
			CONSTRAINT UC_reactant UNIQUE (reactant_reaction, reactant_compound)
		);
		"""

		self.execute(sql)

	def create_table_pose(self):
		"""Create the pose table"""
		logger.debug('HIPPO.Database.create_table_pose()')

		sql = """CREATE TABLE pose(
			pose_id INTEGER PRIMARY KEY,
			pose_inchikey TEXT,
			pose_alias TEXT,
			-- pose_longname TEXT,
			pose_smiles TEXT,
			pose_reference INTEGER,
			pose_path TEXT,
			pose_compound INTEGER,
			pose_target INTEGER,
			pose_mol BLOB,
			pose_fingerprint BLOB,
			pose_energy_score REAL,
			pose_distance_score REAL,
			pose_metadata TEXT,
			FOREIGN KEY (pose_compound) REFERENCES compound(compound_id),
			CONSTRAINT UC_pose_alias UNIQUE (pose_alias)
			CONSTRAINT UC_pose_path UNIQUE (pose_path)
		);
		"""

		### snippet to convert python metadata dictionary with JSON
		# json.dumps(variables).encode('utf-8')
		# json.loads(s.decode('utf-8'))

		self.execute(sql)

	def create_table_tag(self):
		"""Create the tag table"""
		logger.debug('HIPPO.Database.create_table_tag()')

		sql = """CREATE TABLE tag(
			tag_name TEXT,
			tag_compound INTEGER,
			tag_pose INTEGER,
			FOREIGN KEY (tag_compound) REFERENCES compound(compound_id),
			FOREIGN KEY (tag_pose) REFERENCES pose(pose_id),
			CONSTRAINT UC_tag_compound UNIQUE (tag_name, tag_compound)
			CONSTRAINT UC_tag_pose UNIQUE (tag_name, tag_pose)
		);
		"""

		self.execute(sql)

	def create_table_quote(self):
		"""Create the quote table"""
		logger.debug('HIPPO.Database.create_table_quote()')

		sql = """CREATE TABLE quote(
			quote_id INTEGER PRIMARY KEY,
			quote_smiles TEXT,
			quote_amount REAL,
			quote_supplier TEXT,
			quote_catalogue TEXT,
			quote_entry TEXT,
			quote_lead_time INTEGER,
			quote_price REAL,
			quote_currency TEXT,
			quote_purity REAL,
			quote_date TEXT,
			quote_compound INTEGER,
			FOREIGN KEY (quote_compound) REFERENCES compound(compound_id),
			CONSTRAINT UC_quote UNIQUE (quote_amount, quote_supplier, quote_catalogue, quote_entry)
		);
		"""

		self.execute(sql)

		# if logging.root.level >= logging.DEBUG:
		# 	self.execute(f'PRAGMA table_info(quote);')
		# 	pprint(self.cursor.fetchall())

	def create_table_target(self):
		"""Create the target table"""
		logger.debug('HIPPO.Database.create_table_target()')
		sql = """CREATE TABLE target(
			target_id INTEGER PRIMARY KEY,
			target_name TEXT,
			target_metadata TEXT,
			CONSTRAINT UC_target UNIQUE (target_name)
		);
		"""

		self.execute(sql)

	def create_table_feature(self):
		"""Create the feature table"""
		logger.debug('HIPPO.Database.create_table_feature()')
		sql = """CREATE TABLE feature(
			feature_id INTEGER PRIMARY KEY,
			feature_family TEXT,
			feature_target INTEGER,
			feature_chain_name TEXT,
			feature_residue_name TEXT,
			feature_residue_number INTEGER,
			feature_atom_names TEXT,
			CONSTRAINT UC_feature UNIQUE (feature_family, feature_target, feature_chain_name, feature_residue_number, feature_residue_name, feature_atom_names)
		);
		"""

		self.execute(sql)

	def create_table_pattern_bfp(self):
		"""Create the pattern_bfp table"""
		logger.debug('HIPPO.Database.create_table_pattern_bfp()')

		sql = "CREATE VIRTUAL TABLE compound_pattern_bfp USING rdtree(compound_id, fp bits(2048))"
 
		self.execute(sql)

	def create_table_morgan_bfp(self):
		"""Create the morgan_bfp table"""
		logger.debug('HIPPO.Database.create_table_morgan_bfp()')
 
		sql = "CREATE VIRTUAL TABLE compound_morgan_bfp USING rdtree(compound_id, fp bits(2048))"
 
		self.execute(sql)

	### INSERTION

	def insert_compound(self, 
		*,
		smiles: str, 
		base: Compound | int | None = None, 
		alias: str | None = None,
		tags: None | list = None, 
		warn_duplicate=True,
		commit=True,
		metadata=None,
		inchikey: str = None, 
	) -> int:
		"""Insert a compound"""

		# process the base
		assert isinstance(base, int) or isinstance(base, Compound) or base is None, f'incompatible base={base}'

		if base and not isinstance(base,int):
			base = base.id

		# logger.debug(f'{base=}')

		# generate the inchikey name
		inchikey = inchikey or inchikey_from_smiles(smiles)

		# logger.debug(f'{smiles} --> {name}')

		sql = """
		INSERT INTO compound(compound_inchikey, compound_smiles, compound_base, compound_mol, compound_pattern_bfp, compound_morgan_bfp, compound_alias)
		VALUES(?1, ?2, ?3, mol_from_smiles(?2), mol_pattern_bfp(mol_from_smiles(?2), 2048), mol_morgan_bfp(mol_from_smiles(?2), 2, 2048), ?4)
		"""

		try:
			self.execute(sql, (inchikey, smiles, base, alias))

		except sqlite3.IntegrityError as e:
			if 'UNIQUE constraint failed: compound.compound_inchikey' in str(e):
				if warn_duplicate:
					logger.warning(f"Skipping compound with existing inchikey \"{inchikey}\"")
			elif 'UNIQUE constraint failed: compound.compound_smiles' in str(e):
				if warn_duplicate:
					logger.warning(f"Skipping compound with existing smiles \"{smiles}\"")
			elif 'UNIQUE constraint failed: compound.compound_pattern_bfp' in str(e):
				if warn_duplicate:
					logger.warning(f"Skipping compound with existing pattern binary fingerprint \"{smiles}\"")
			elif 'UNIQUE constraint failed: compound.compound_morgan_bfp' in str(e):
				if warn_duplicate:
					logger.warning(f"Skipping compound with existing morgan binary fingerprint \"{smiles}\"")
			else:
				logger.exception(e)

			# if tags:
				# logger.error('Unused tag data')
			return None

		except Exception as e:
			logger.exception(e)

		compound_id = self.cursor.lastrowid
		if commit:
			self.commit()

		### register the tags

		if tags:
			for tag in tags:
				self.insert_tag(name=tag, compound=compound_id, commit=commit)

		### register the binary fingerprints

		result = self.insert_compound_pattern_bfp(compound_id, commit=commit)

		if not result:
			logger.error('Could not insert compound pattern bfp')

		# ### register the morgan fingerprint

		# result = self.insert_compound_morgan_bfp(compound_id)

		# if not result:
		# 	logger.error('Could not insert compound morgan bfp')

		if metadata:
			self.insert_metadata(table='compound', id=compound_id, payload=metadata, commit=commit)

		return compound_id

	def insert_compound_pattern_bfp(self, compound_id, commit=True):
		"""Insert a compound_pattern_bfp"""

		sql = """
		INSERT INTO compound_pattern_bfp(compound_id, fp)
		VALUES(?1, ?2)
		"""
		
		bfp, = self.select_where('compound_pattern_bfp', 'compound', 'id', compound_id)

		try:
			self.execute(sql, (compound_id, bfp))

		except Exception as e:
			logger.exception(e)

		bfp_id = self.cursor.lastrowid
		if commit:
			self.commit()

		return bfp_id

	def insert_compound_morgan_bfp(self, compound_id):
		"""Insert a compound_morgan_bfp"""


		sql = """
		INSERT INTO compound_morgan_bfp(compound_id, fp)
		VALUES(?1, ?2)
		"""

		bfp, = self.select_where('compound_morgan_bfp', 'compound', 'id', compound_id)

		try:
			self.execute(sql, (compound_id, bfp))

		except Exception as e:
			logger.exception(e)

		mfp_id = self.cursor.lastrowid
		self.commit()

		return mfp_id

	def insert_pose(self,
		*,
		compound: Compound | int,
		target: int | str,
		path: str,
		inchikey: str | None = None,
		alias: str | None = None,
		reference: int | Pose | None = None,
		tags: None | list = None,
		energy_score: float | None = None,
		distance_score: float | None = None,
		metadata: None | dict = None,
		commit: bool = True,
		longname: str = None,
		warn_duplicate: bool = True,
		resolve_path: bool = True,
	):
		"""Insert a pose"""

		if isinstance(compound, Compound):
			compound = compound.id

		if isinstance(reference, Pose):
			reference = reference.id

		# assert isinstance(compound, Compound), f'incompatible {compound}'
		# assert isinstance(reference, Pose) or reference is None, f'incompatible reference={reference}'

		if isinstance(target, str):
			target = self.get_target_id(name=target)
			if not target:
				raise ValueError(f'No such {target=}')
		target_name = self.get_target_name(id=target)
		
		if resolve_path:
			try:
				path = Path(path)
				path = path.resolve(strict=True)
				path = str(path)

			except FileNotFoundError as e:
				logger.error(f'Path cannot be resolved: {mcol.file}{path}')
				return None

		sql = """
		INSERT INTO pose(pose_inchikey, pose_alias, pose_smiles, pose_compound, pose_target, pose_path, pose_reference, pose_energy_score, pose_distance_score)
		VALUES(?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9)
		"""

		try:
			self.execute(sql, (inchikey, alias, None, compound, target, path, reference, energy_score, distance_score))

		except sqlite3.IntegrityError as e:
			if 'UNIQUE constraint failed: pose.pose_path' in str(e):
				if warn_duplicate:
					logger.warning(f"Could not insert pose with duplicate path \"{path}\"")
			elif 'UNIQUE constraint failed: pose.pose_alias' in str(e):
				if warn_duplicate:
					logger.warning(f"Could not insert pose with duplicate alias \"{alias}\"")
			else:
				logger.exception(e)
			return None

		except Exception as e:
			logger.exception(e)
			raise
			
		pose_id = self.cursor.lastrowid
		if commit:
			self.commit()

		if tags:
			for tag in tags:
				self.insert_tag(name=tag, pose=pose_id, commit=commit)

		if metadata:
			self.insert_metadata(table='pose', id=pose_id, payload=metadata, commit=commit)

		return pose_id

	def insert_tag(self, 
		*,
		name, 
		compound: int = None, 
		pose: int = None,
		commit: bool = True,
	):
		"""Insert a tag"""

		assert compound or pose

		sql = """
		INSERT INTO tag(tag_name, tag_compound, tag_pose)
		VALUES(?1, ?2, ?3)
		"""

		try:
			self.execute(sql, (name, compound, pose))

		except sqlite3.IntegrityError as e:
			return None

		except Exception as e:
			logger.exception(e)		

		tag_id = self.cursor.lastrowid
		if commit:
			self.commit()
		return tag_id

	def insert_inspiration(self, 
		*,
		original: Pose | int, 
		derivative: Pose | int,
		warn_duplicate: bool = True,
		commit: bool = True,
	) -> int:
		"""Insert an inspiration"""

		if isinstance(original, Pose):
			original = original.id
		if isinstance(derivative, Pose):
			derivative = derivative.id

		assert isinstance(original, int), 'Must pass an integer ID or Pose object (original)'
		assert isinstance(derivative, int), 'Must pass an integer ID or Pose object (derivative)'

		sql = """
		INSERT INTO inspiration(inspiration_original, inspiration_derivative)
		VALUES(?1, ?2)
		"""

		try:
			self.execute(sql, (original, derivative))

		except sqlite3.IntegrityError as e:
			if warn_duplicate:
				logger.warning(f"Skipping existing inspiration: {original} {derivative}")
			return None

		except Exception as e:
			logger.exception(e)

		inspiration_id = self.cursor.lastrowid

		if commit:
			self.commit()
		return inspiration_id

	def insert_reaction(self,
		*,
		type: str,
		product: Compound | int,
		product_yield: float = 1.0,
		commit: bool = True,
	) -> int:
		"""Insert a reaction"""

		if isinstance(product, Compound):
			product = product.id

		# assert isinstance(product, Compound), f'incompatible {product=}'
		assert isinstance(type, str), f'incompatible {type=}'

		sql = """
		INSERT INTO reaction(reaction_type, reaction_product, reaction_product_yield)
		VALUES(?1, ?2, ?3)
		"""

		try:
			self.execute(sql, (type, product, product_yield))

		except Exception as e:
			logger.exception(e)
			
		reaction_id = self.cursor.lastrowid
		if commit:
			self.commit()
		return reaction_id

	def insert_reactant(self,
		*,
		compound: Compound,
		reaction: Reaction,
		amount: float = 1.0,
		commit: bool = True,
	) -> int:
		"""Insert a reactant"""

		if isinstance(reaction, int):
			reaction = self.get_reaction(id=reaction)

		if isinstance(compound, int):
			compound = self.get_compound(id=compound)

		assert isinstance(compound, Compound), f'incompatible {compound=}'
		assert isinstance(reaction, Reaction), f'incompatible {reaction=}'

		sql = """
		INSERT INTO reactant(reactant_amount, reactant_reaction, reactant_compound)
		VALUES(?1, ?2, ?3)
		"""

		try:
			self.execute(sql, (amount, reaction.id, compound.id))

		except sqlite3.IntegrityError as e:
			logger.warning(f"Skipping existing reactant: {reaction=} {compound=}")

		except Exception as e:
			logger.exception(e)
			
		reactant_id = self.cursor.lastrowid
		if commit:
			self.commit()
		return reactant_id

	def insert_quote(self,
		*,
		compound: Compound | int,
		supplier: str,
		catalogue: str | None,
		entry: str | None,
		amount: float,
		price: float,
		currency: str | None,
		purity: float | None,
		lead_time: int,
		smiles: str | None = None,
		commit: bool = True,
	) -> int | None:
		"""Insert a quote"""

		if not isinstance(compound, int):
			assert isinstance(compound, Compound), f'incompatible {compound=}'
			compound = compound.id
			
		assert currency in ['GBP', 'EUR', 'USD', None], f'incompatible {currency=}'
		
		assert supplier in ['MCule', 'Enamine', 'Stock'], f'incompatible {supplier=}'

		smiles = smiles or ""
		
		sql = """
		INSERT or REPLACE INTO quote(
			quote_smiles,
			quote_amount,
			quote_supplier,
			quote_catalogue,
			quote_entry,
			quote_lead_time,
			quote_price,
			quote_currency,
			quote_purity,
			quote_compound,
			quote_date
		)
		VALUES(?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, date());
		"""

		try:
			self.execute(sql, (smiles, amount, supplier, catalogue, entry, lead_time, price, currency, purity, compound))

		except sqlite3.InterfaceError as e:
			logger.error(e)
			logger.debug((smiles, amount, supplier, catalogue, entry, lead_time, price, currency, purity, compound))
			raise

		except Exception as e:
			logger.exception(e)
			return None
			
		quote_id = self.cursor.lastrowid
		if commit:
			self.commit()
		return quote_id

	def insert_target(self,
		*,
		name: str
	) -> int:
		"""Insert a target"""

		sql = """
		INSERT INTO target(target_name)
		VALUES(?1)
		"""

		try:
			self.execute(sql, (name, ))

		except sqlite3.IntegrityError as e:
			logger.warning(f"Skipping existing target with {name=}")
			return None

		except Exception as e:
			logger.exception(e)
			
		target_id = self.cursor.lastrowid
		self.commit()
		return target_id

	def insert_feature(self,
		*,
		family: str,
		target: int,
		chain_name: str,
		residue_name: str,
		residue_number: int,
		atom_names: list,
	) -> int:
		"""Insert a feature"""

		assert len(chain_name) == 1
		assert len(residue_name) <= 4
		for a in atom_names:
			assert len(a) <= 4

		if isinstance(target, str):
			target = self.get_target_id(name=target)
		assert isinstance(target, int)

		sql = """
		INSERT INTO feature(feature_family, feature_target, feature_chain_name, feature_residue_name, feature_residue_number, feature_atom_names)
		VALUES(?1, ?2, ?3, ?4, ?5, ?6)
		"""

		atom_names = ' '.join(atom_names)

		try:
			self.execute(sql, (family, target, chain_name, residue_name, residue_number, atom_names))

		except sqlite3.IntegrityError as e:
			return None

		except Exception as e:
			logger.exception(e)
			
		target_id = self.cursor.lastrowid
		self.commit()
		return target_id

	def insert_metadata(self,
		*,
		table: str,
		id: int,
		payload: dict,
		commit: bool = True,
	) -> None:
		"""Insert metadata"""

		payload = json.dumps(payload)

		self.update(table=table, id=id, key=f'{table}_metadata', value=payload, commit=commit)

	### SELECTION

	def select(self, query, table, multiple=False):
		"""Select entries"""

		sql = f'SELECT {query} FROM {table}'

		try:
			self.execute(sql)
		except sqlite3.OperationalError as e:
			logger.var('sql',sql)
			raise

		if multiple:
			result = self.cursor.fetchall()
		else:
			result = self.cursor.fetchone()

		return result

	def select_where(self, query, table, key, value=None, multiple=False, none='error', sort=None):
		"""Select entries where key==value"""

		if isinstance(value, str):
			value = f"'{value}'"

		if value is not None:
			where_str = f'{table}_{key}={value}'
		else:
			where_str = key

		if sort:
			sql = f'SELECT {query} FROM {table} WHERE {where_str} ORDER BY {sort}'
		else:
			sql = f'SELECT {query} FROM {table} WHERE {where_str}'

		try:
			self.execute(sql)
		except sqlite3.OperationalError as e:
			logger.var('sql',sql)
			raise

		if multiple:
			result = self.cursor.fetchall()
		else:
			result = self.cursor.fetchone()

		if not result and none == 'error':
			logger.error(f'No entry in {table} with {where_str}')
			return None

		return result

	def select_id_where(self, table, key, value, multiple=False, none='error'):
		"""Select ID's where key==value"""
		return self.select_where(query=f'{table}_id', table=table, key=key, value=value, multiple=multiple, none=none)

	def select_all_where(self, table, key, value, multiple=False, none='error'):
		"""Select * where key==value"""
		return self.select_where(query='*', table=table, key=key, value=value, multiple=multiple, none=none)

	### DELETION

	def delete_where(self, table, key, value=None):
		"""Delete where key==value"""

		if value is not None:

			if isinstance(value, str):
				value = f"'{value}'"

			sql = f'DELETE FROM {table} WHERE {table}_{key}={value}'

		else:
			
			sql = f'DELETE FROM {table} WHERE {key}'

		try:
			result = self.execute(sql)
			self.commit()

		except sqlite3.OperationalError as e:
			logger.var('sql',sql)
			raise

		return None

	def delete_tag(self, tag):
		self.delete_where(table='tag', key='name', value=tag)

	### UPDATE

	def update(self, *, 
		table: str, 
		id: int, 
		key: str, 
		value, 
		commit: bool = True,
	):
		"""Update a database entry where key==value"""

		# sql = f"""
		# UPDATE ?1
		# SET ?2 = ?3
		# WHERE ?4 = ?5;
		# """
		
		sql = f"""
		UPDATE {table}
		SET {key} = ?
		WHERE {table}_id = {id};
		"""
		
		try:
			# self.execute(sql, (table, key, value, f'{table}_id', id))
			self.execute(sql, (value, ))
		except sqlite3.OperationalError as e:
			logger.var('sql',sql)
			raise

		id = self.cursor.lastrowid
		if commit:
			self.commit()
		return id

	### GETTERS

	def get_compound(self,
		*,
		table: str = 'compound',
		id: int | None = None,
		inchikey: str | None = None,
		alias: str | None = None,
		smiles: str | None = None,
	) -> Compound:
		"""Get a compound"""
		
		if id is None:
			id = self.get_compound_id(inchikey=inchikey, smiles=smiles, alias=alias)

		if not id:
			logger.error(f'Invalid {id=}')
			return None

		query = 'compound_id, compound_inchikey, compound_alias, compound_smiles, compound_base'
		entry = self.select_where(query=query, table=table, key='id', value=id)
		compound = Compound(self._animal, self, *entry, metadata=None, mol=None)
		return compound

	def get_compound_id(self, 
		*,
		table: str = 'compound',
		inchikey: str | None = None, 
		alias: str | None = None, 
		smiles: str | None = None, 
	) -> int:
		"""Get a compound ID"""

		if inchikey:
			entry = self.select_id_where(table=table, key='inchikey', value=inchikey)

		elif alias:
			entry = self.select_id_where(table=table, key='alias', value=alias)

		elif smiles:
			entry = self.select_id_where(table=table, key='smiles', value=smiles)

		else:
			raise NotImplementedError

		if entry:
			return entry[0]

		return None

	def get_compound_computed_property(self, prop, compound_id):
		"""Use chemicalite to calculate a property from the stored binary molecule"""
		function = CHEMICALITE_COMPOUND_PROPERTY_MAP[prop]
		val, = self.select_where(query=f'{function}(compound_mol)', table='compound', key='id', value=compound_id, multiple=False)
		return val

	def get_pose(self,
		*,
		table: str = 'pose',
		id: int | None = None,
		inchikey: str = None,
		alias: str = None,
	) -> Pose:
		"""Get a pose"""
		
		if id is None:
			id = self.get_pose_id(inchikey=inchikey, alias=alias)

		if not id:
			logger.error(f'Invalid {id=}')
			return None

		query = 'pose_id, pose_inchikey, pose_alias, pose_smiles, pose_reference, pose_path, pose_compound, pose_target, pose_mol, pose_fingerprint, pose_energy_score, pose_distance_score'
		entry = self.select_where(query=query, table=table, key='id', value=id)
		pose = Pose(self, *entry)
		return pose

	def get_pose_id(self, 
		*,
		table: str = 'pose',
		inchikey: str | None = None, 
		alias: str | None = None, 
	) -> int:
		"""Get a pose ID"""

		if inchikey:
			# inchikey might not be unique
			entries = self.select_id_where(table=table, key='inchikey', value=inchikey, multiple=True)
			if len(entries) != 1:
				logger.warning(f'Multiple poses with {inchikey=}')
				return entries
			else:
				entry = entries[0]

		elif alias:
			entry = self.select_id_where(table=table, key='alias', value=alias)
		
		else:
			raise NotImplementedError
			
		if entry:
			return entry[0]

		return None

	def get_reaction(self,
		*,
		table: str = 'reaction',
		id: int | None = None,
		none: str | None = None,
	) -> Reaction:
		"""Get a reaction"""
		
		if not id:
			logger.error(f'Invalid {id=}')
			return None

		query = 'reaction_id, reaction_type, reaction_product, reaction_product_yield'
		entry = self.select_where(query=query, table=table, key='id', value=id, none=none)
		reaction = Reaction(self, *entry)
		return reaction

	def get_quote(self,
		*,
		table: str = 'quote',
		id: int | None = None,
		none: str | None = None,
	) -> list[dict]:
		"""Get a quote"""
		
		query = 'quote_compound, quote_supplier, quote_catalogue, quote_entry, quote_amount, quote_price, quote_currency, quote_lead_time, quote_purity, quote_date, quote_smiles, quote_id '
		entry = self.select_where(query=query, table=table, key='id', value=id, none=none)

		return Quote(
			db=self,
			id=entry[11],
			compound=entry[0],
			supplier=entry[1],
			catalogue=entry[2],
			entry=entry[3],
			amount=entry[4],
			price=entry[5],
			currency=entry[6],
			lead_time=entry[7],
			purity=entry[8],
			date=entry[9],
			smiles=entry[10],
		)

	def get_metadata(self, 
		*,
		table: str, 
		id: int
	) -> dict:
		"""Get metadata"""

		payload, = self.select_where(query=f'{table}_metadata', table=table, key=f'id', value=id)

		if payload:
			payload = json.loads(payload)
			
		else:
			payload = dict()

		metadata = MetaData(payload)

		metadata._db = self
		metadata._id = id
		metadata._table = table

		return metadata

	def get_target(self, 
		*,
		id=int
	) -> Target:
		"""Get target"""
		return Target(db=self,id=id,name=self.get_target_name(id=id))

	def get_target_name(self,
		*,
		id:int,
	) -> str:
		"""Get target name"""

		table = 'target'
		payload, = self.select_where(query=f'{table}_name', table=table, key='id', value=id)
		return payload

	def get_target_id(self, 
		*,
		name: str, 
	) -> int:
		"""Get target ID"""
			
		table = 'target'
		entry = self.select_id_where(table=table, key='name', value=name)

		if entry:
			return entry[0]

		return None

	def get_feature(self, 
		*,
		id=int
	) -> Feature:
		"""Get feature"""
		entry = self.select_all_where(table='feature', key='id', value=id)
		return Feature(*entry)

	def get_possible_reaction_ids(self, *, compound_ids):

		compound_ids_str = str(tuple(compound_ids)).replace(',)',')')

		result = self.execute(f'''
			WITH possible_reactants AS 
			(
				SELECT reactant_reaction, CASE WHEN reactant_compound IN {compound_ids_str} THEN reactant_compound END AS [possible_reactant] FROM reactant
			)

			, possible_reactions AS (
				SELECT reactant_reaction, COUNT(CASE WHEN possible_reactant IS NULL THEN 1 END) AS [count_null] FROM possible_reactants
				GROUP BY reactant_reaction
			)
			
			SELECT reactant_reaction FROM possible_reactions
			WHERE count_null = 0
		''').fetchall()

		return [q for q, in result]

	def get_unsolved_reaction_tree(self, *, product_ids, debug=False):

		from .cset import CompoundSet
		from .rset import ReactionSet

		all_reactants = set()
		all_reactions = set()

		for i in range(300):

			if debug: 
				logger.var('recursive depth',i+1)
			
			if debug: 
				logger.var('#products', len(product_ids))

			product_ids_str = str(tuple(product_ids)).replace(',)',')')

			reaction_ids = self.select_where(table='reaction', query='DISTINCT reaction_id', key=f'reaction_product in {product_ids_str}', multiple=True, none='quiet')
			
			reaction_ids = [q for q, in reaction_ids]

			if not reaction_ids:
				break

			for reaction_id in reaction_ids:
				all_reactions.add(reaction_id)

			if debug: 
				logger.var('#reactions', len(reaction_ids))

			reaction_ids_str = str(tuple(reaction_ids)).replace(',)',')')

			reactant_ids = self.select_where(table='reactant', query='DISTINCT reactant_compound', key=f'reactant_reaction in {reaction_ids_str}', multiple=True)
			
			if debug: 
				logger.var('#reactants', len(reactant_ids))
			
			reactant_ids = [q for q, in reactant_ids]

			if not reactant_ids:
				break

			for reactant_id in reactant_ids:
				all_reactants.add(reactant_id)

			product_ids = reactant_ids

		# all intermediates
		ids = self.execute('SELECT DISTINCT reaction_product FROM reaction INNER JOIN reactant ON reaction.reaction_product = reactant.reactant_compound').fetchall()
		ids = [q for q, in ids]
		intermediates = CompoundSet(self, ids)

		# remove intermediates
		cset = CompoundSet(self, all_reactants)
		all_reactants = cset - intermediates

		# reactions
		all_reactions = ReactionSet(self, all_reactions)

		if debug:
			logger.var('#all_reactants', len(all_reactants))
			logger.var('#all_reactions', len(all_reactions))

		return all_reactants, all_reactions

	def get_possible_reaction_product_ids(self, *, reaction_ids):
		reaction_ids_str = str(tuple(reaction_ids)).replace(',)',')')
		return [q for q, in self.select_where(query='DISTINCT reaction_product', table='reaction', key=f'reaction_id IN {reaction_ids_str}', multiple=True)] 

	### COMPOUND QUERY

	def query_substructure(self, query, fast=True, none='error'):
		"""Search compounds by substructure"""

		# smiles
		if isinstance(query, str):

			if fast:
				sql = f"SELECT compound.compound_id FROM compound, compound_pattern_bfp AS bfp WHERE compound.compound_id = bfp.compound_id AND mol_is_substruct(compound.compound_mol, mol_from_smiles(?))"
			else:
				sql = f"SELECT compound_id FROM compound WHERE mol_is_substruct(compound_mol, mol_from_smiles(?))"

		else:

			raise NotImplementedError

		try:
			self.execute(sql, (query, ))
		except sqlite3.OperationalError as e:
			logger.var('sql',sql)
			raise

		result = self.cursor.fetchall()

		if not result and none == 'error':
			logger.error(f'No compounds with substructure {query}')
			return None

		compounds = []
		for id, in result:
			compounds.append(self.get_compound(id=id))

		return compounds

	def query_similarity(self, query, threshold, return_similarity=False, none='error'):
		"""Search compounds by similarity"""

		from .cset import CompoundSet

		# smiles
		if isinstance(query, str):

			if return_similarity:
				# sql = f"SELECT compound_id, bfp_tanimoto(mol_morgan_bfp(mol_from_smiles(?1), 2, 1024), mol_morgan_bfp(compound.compound_mol, 2, 1024)) as t FROM compound JOIN compound_morgan_bfp AS mfp USING(compound_id) WHERE mfp.compound_id match rdtree_tanimoto(mol_morgan_bfp(mol_from_smiles(?1), 2, 1024), ?2) ORDER BY t DESC "
				sql = f"SELECT compound_id, bfp_tanimoto(mol_pattern_bfp(mol_from_smiles(?1), 2048), mol_pattern_bfp(compound.compound_mol, 2048)) as t FROM compound JOIN compound_pattern_bfp AS mfp USING(compound_id) WHERE mfp.compound_id match rdtree_tanimoto(mol_pattern_bfp(mol_from_smiles(?1), 2048), ?2) ORDER BY t DESC "
			else:
				# sql = f"SELECT compound_id FROM compound_morgan_bfp AS mfp WHERE mfp.compound_id match rdtree_tanimoto(mol_morgan_bfp(mol_from_smiles(?1), 2, 1024), ?2) "
				sql = f"SELECT compound_id FROM compound_pattern_bfp AS mfp WHERE mfp.compound_id match rdtree_tanimoto(mol_pattern_bfp(mol_from_smiles(?1), 2048), ?2) "

		else:
			raise NotImplementedError

		try:
			self.execute(sql, (query, threshold))
		except sqlite3.OperationalError as e:
			logger.var('sql',sql)
			raise

		result = self.cursor.fetchall()

		if not result and none == 'error':
			logger.error(f'No compounds with similarity >= {threshold} to {query}')
			return None

		if return_similarity:
			ids, similarities = zip(*result)
			cset = CompoundSet(self, ids)
			return cset, similarities

		ids = [r for r, in result]
		cset = CompoundSet(self, ids)

		return cset

	def query_exact(self, query):
		"""Search for exact match compounds"""
		return self.query_similarity(query, 0.989, return_similarity=False)

	### COUNTING

	def count(self, table):
		"""Count all entries in a table"""

		sql = f"""SELECT COUNT(1) FROM {table}; """
		self.execute(sql)
		return self.cursor.fetchone()[0]

	def count_where(self, table, key, value=None):
		"""Count all entries in a table where key==value"""
		if value is not None:
			where_str = f'{table}_{key} is {value}'
		else:
			where_str = key

		sql = f"""SELECT COUNT(1) FROM {table} WHERE {where_str};"""
		self.execute(sql)
		return self.cursor.fetchone()[0]

	### ID SELECTION

	def min_id(self, table):
		"""Return the smallest ID in the given table"""
		id, = self.select(table=table, query=f'MIN({table}_id)')
		return id

	def max_id(self, table):
		"""Return the largest ID in the given table"""
		id, = self.select(table=table, query=f'MAX({table}_id)')
		return id

	def slice_ids(self, *, table, start, stop, step=1):
		"""Return a list of indices in the given slice"""

		min_id = self.min_id(table)
		max_id = self.max_id(table)

		start = start or min_id
		stop = stop or max_id + 1
		step = step or 1

		if not (start >= 0 and start <= max_id):
			raise IndexError(f'Slice {start=} outside of DB {table}_id range ({min_id}, {max_id})')

		if not (stop >= 0 and stop <= max_id + 1):
			raise IndexError(f'Slice {stop=} outside of DB {table}_id range ({min_id}, {max_id})')

		if step != 1:
			raise NotImplementedError(f'Slice {step=} not supported')

		ids = self.select_where(table=table, query=f'{table}_id', key=f'{table}_id >= {start} AND {table}_id < {stop}', multiple=True)
		ids = [q for q, in ids]

		return ids

	### PRUNING

	def prune_reactions(self, compound, reactions):
		"""Remove duplicate reactions for a given compound"""

		pruned = []
		del_list = []

		for i,reaction in enumerate(reactions):

			matches = [r for r in pruned if r == reaction]

			if not matches:
				pruned.append(reaction)

			else:
				del_list.append(reaction)

		for reaction in del_list:
			logger.warning(f'Deleted duplicate {reaction=}')
			self.delete_where('reaction', 'id', reaction.id)

		return pruned

	def remove_metadata_list_item(self, *, table, key, value, remove_empty=True):
		"""Remove a specific item from list-like values associated with a given key from all metadata entries in a given table"""

		# get id's with specific metadata key and value
		value_str = json.dumps(value)
		result = self.select_where(query=f'{table}_id, {table}_metadata', table=table, key=f'{table}_metadata LIKE \'%"export": [%{value_str}%]%\'', multiple=True, none='quiet')

		# loop over all matches
		for id, metadata_str in result:

			# read the metadata
			metadata = json.loads(metadata_str)
			
			# modify the metadata
			index = metadata[key].index(value)
			metadata[key].pop(index)

			if remove_empty and not metadata[key]:
				del metadata[key]
			
			# update the database
			metadata_str = json.dumps(metadata)
			self.update(table=table, id=id, key=f'{table}_metadata', value=metadata_str, commit=False)

		# only runs if non-zero matches
		else:
			# commit the changes
			self.commit()

	### PRINTING

	def list_tables(self):
		"""List all tables in the database"""
		self.execute("SELECT name FROM sqlite_master WHERE type='table';")
		pprint(self.cursor.fetchall())

	def print_table(self, name):
		"""Print a table's entries"""
		self.execute(f"SELECT * FROM {name}")
		pprint(self.cursor.fetchall())

	### DUNDERS

CHEMICALITE_COMPOUND_PROPERTY_MAP = {
	'num_heavy_atoms':'mol_num_hvyatms',
	'formula':'mol_formula',
	'num_rings':'mol_num_rings',
}
