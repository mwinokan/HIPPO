
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

	def __init__(self, path: Path) -> None:

		assert isinstance(path, Path)
		
		logger.debug('hippo3.Database.__init__()')

		self._path = path
		self._connection = None
		self._cursor = None

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
		return self._path

	@property
	def connection(self):
		if not self._connection:
			self.connect()
		return self._connection
	
	@property
	def cursor(self):
		if not self._cursor:
			self._cursor = self.connection.cursor()
		return self._cursor

	### PUBLIC METHODS / API CALLS

	def close(self) -> None:
		logger.debug('hippo3.Database.close()')
		if self.connection:
			self.connection.close()
		mout.success(f'Closed connection to {mcol.file}{self.path}')

	### GENERAL SQL

	def connect(self):
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
		try:
			if payload:
				return self.cursor.execute(sql, payload)
			else:
				return self.cursor.execute(sql)
		except Error as e:
			# logger.exception(e)
			raise

	def commit(self):
		# logger.debug('db.commit...')
		self.connection.commit()
		# raise Exception

	### CREATE TABLES

	def create_blank_db(self):

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
		logger.debug('HIPPO.Database.create_table_compound()')

		sql = """CREATE TABLE compound(
			compound_id INTEGER PRIMARY KEY,
			compound_name TEXT,
			compound_smiles TEXT,
			compound_base INTEGER,
			compound_mol MOL,
			compound_pattern_bfp bits(2048),
			compound_morgan_bfp bits(2048),
			compound_metadata TEXT,
			FOREIGN KEY (compound_base) REFERENCES compound(compound_id),
			CONSTRAINT UC_compound_name UNIQUE (compound_name)
			CONSTRAINT UC_compound_smiles UNIQUE (compound_smiles)
		);
		"""
			# CONSTRAINT UC_compound_morgan_bfp UNIQUE (compound_morgan_bfp)
			# CONSTRAINT UC_compound_pattern_bfp UNIQUE (compound_pattern_bfp)

		self.execute(sql)

	def create_table_inspiration(self):
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
		logger.debug('HIPPO.Database.create_table_pose()')

		sql = """CREATE TABLE pose(
			pose_id INTEGER PRIMARY KEY,
			pose_name TEXT,
			pose_longname TEXT,
			pose_smiles TEXT,
			pose_reference INTEGER,
			pose_path TEXT,
			pose_compound INTEGER,
			pose_target INTEGER,
			pose_mol BLOB,
			pose_fingerprint BLOB,
			pose_metadata TEXT,
			FOREIGN KEY (pose_compound) REFERENCES compound(compound_id),
			CONSTRAINT UC_pose_longname UNIQUE (pose_longname)
			CONSTRAINT UC_pose_path UNIQUE (pose_path)
		);
		"""

		### snippet to convert python metadata dictionary with JSON
		# json.dumps(variables).encode('utf-8')
		# json.loads(s.decode('utf-8'))

		self.execute(sql)

	def create_table_tag(self):
		logger.debug('HIPPO.Database.create_table_tag()')

		sql = """CREATE TABLE tag(
			tag_name TEXT,
			tag_compound INTEGER,
			tag_pose INTEGER,
			FOREIGN KEY (tag_compound) REFERENCES compound(compound_id),
			FOREIGN KEY (tag_pose) REFERENCES pose(pose_id),
			CONSTRAINT UC_tag UNIQUE (tag_name, tag_compound, tag_pose)
		);
		"""

		self.execute(sql)

	def create_table_quote(self):
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
		logger.debug('HIPPO.Database.create_table_feature()')
		sql = """CREATE TABLE feature(
			feature_id INTEGER PRIMARY KEY,
			feature_family TEXT,
			feature_target INTEGER,
			feature_chain_name TEXT,
			feature_residue_name TEXT,
			feature_residue_number INTEGER,
			feature_atom_names TEXT,
			CONSTRAINT UC_feature UNIQUE (feature_family, feature_target, feature_chain_name, feature_residue_number, feature_atom_names)
		);
		"""

		self.execute(sql)

	def create_table_pattern_bfp(self):
		logger.debug('HIPPO.Database.create_table_pattern_bfp()')

		sql = "CREATE VIRTUAL TABLE compound_pattern_bfp USING rdtree(compound_id, fp bits(2048))"
 
		self.execute(sql)

	def create_table_morgan_bfp(self):
		logger.debug('HIPPO.Database.create_table_morgan_bfp()')
 
		sql = "CREATE VIRTUAL TABLE compound_morgan_bfp USING rdtree(compound_id, fp bits(2048))"
 
		self.execute(sql)

	### INSERTION

	def insert_compound(self, 
		*,
		smiles: str, 
		base: Compound | int | None = None, 
		tags: None | list = None, 
		warn_duplicate=True,
		commit=True,
		metadata=None,
		inchikey: str = None, 
	) -> int:

		# process the base
		assert isinstance(base, int) or isinstance(base, Compound) or base is None, f'incompatible base={base}'

		if base and not isinstance(base,int):
			base = base.id

		# logger.debug(f'{base=}')

		# generate the inchikey name
		if inchikey:
			name = inchikey
		else:
			name = inchikey_from_smiles(smiles)

		# logger.debug(f'{smiles} --> {name}')

		sql = """
		INSERT INTO compound(compound_name, compound_smiles, compound_base, compound_mol, compound_pattern_bfp, compound_morgan_bfp)
		VALUES(?1, ?2, ?3, mol_from_smiles(?2), mol_pattern_bfp(mol_from_smiles(?2), 2048), mol_morgan_bfp(mol_from_smiles(?2), 2, 2048))
		"""

		try:
			self.execute(sql, (name, smiles, base))

		except sqlite3.IntegrityError as e:
			if 'UNIQUE constraint failed: compound.compound_name' in str(e):
				if warn_duplicate:
					logger.warning(f"Skipping compound with existing name \"{name}\"")
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
		name: str,
		target: int | str,
		path: str,
		reference: int | None = None,
		tags: None | list = None,
		metadata: None | dict = None,
		commit: bool = True,
		longname: str = None,
		warn_duplicate: bool = True,
		resolve_path: bool = True,
	):

		if isinstance(compound, int):
			compound = self.get_compound(id=compound)

		if isinstance(reference, int):
			reference = self.get_pose(id=reference)

		assert isinstance(compound, Compound), f'incompatible {compound}'
		assert name, f'incompatible name={name}'
		assert isinstance(reference, Pose) or reference is None, f'incompatible pose={name}'

		if reference:
			reference = reference.id

		if isinstance(target, str):
			target = self.get_target_id(name=target)
		target_name = self.get_target_name(id=target)
		
		if not longname:
			longname = f'{compound.name} {target_name} {name}'

		if resolve_path:
			try:
				path = Path(path)
				path = path.resolve(strict=True)
				path = str(path)

			except FileNotFoundError as e:
				logger.error(f'Path cannot be resolved: {mcol.file}{path}')
				return None

		sql = """
		INSERT INTO pose(pose_name, pose_longname, pose_smiles, pose_compound, pose_target, pose_path, pose_reference)
		VALUES(?1, ?2, ?3, ?4, ?5, ?6, ?7)
		"""

		try:
			self.execute(sql, (name, longname, compound.smiles, compound.id, target, path, reference))

		except sqlite3.IntegrityError as e:
			if 'UNIQUE constraint failed: pose.pose_longname' in str(e):
				if warn_duplicate:
					logger.warning(f"Skipping pose with existing longname \"{longname}\"")
			elif 'UNIQUE constraint failed: pose.pose_path' in str(e):
				if warn_duplicate:
					logger.warning(f"Skipping pose with existing path \"{longname}\"")
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

		assert compound or pose

		sql = """
		INSERT INTO tag(tag_name, tag_compound, tag_pose)
		VALUES(?1, ?2, ?3)
		"""

		try:
			self.execute(sql, (name, compound, pose))

		except sqlite3.IntegrityError as e:
			if 'UNIQUE constraint failed: compound.compound_name' in str(e):
				logger.warning(f"Skipping compound with existing name \"{name}\"")
			elif 'UNIQUE constraint failed: compound.compound_smiles' in str(e):
				logger.warning(f"Skipping compound with existing smiles \"{smiles}\"")
			else:
				logger.exception(e)

		except Exception as e:
			logger.exception(e)		

		tag_id = self.cursor.lastrowid
		if commit:
			self.commit()
		return tag_id

	def insert_inspiration(self, 
		*,
		original: Pose, 
		derivative: Pose,
		warn_duplicate: bool = True,
		commit: bool = True,
	) -> int:

		assert isinstance(original, Pose), f'incompatible {original=}'
		assert isinstance(derivative, Pose), f'incompatible {derivative=}'

		sql = """
		INSERT INTO inspiration(inspiration_original, inspiration_derivative)
		VALUES(?1, ?2)
		"""

		try:
			self.execute(sql, (original.id, derivative.id))

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
		product: Compound,
		product_yield: float = 1.0,
		commit: bool = True,
	) -> int:

		if isinstance(product, int):
			product = self.get_compound(id=product)

		assert isinstance(product, Compound), f'incompatible {product=}'
		assert isinstance(type, str), f'incompatible {type=}'

		sql = """
		INSERT INTO reaction(reaction_type, reaction_product, reaction_product_yield)
		VALUES(?1, ?2, ?3)
		"""

		try:
			self.execute(sql, (type, product.id, product_yield))

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
		compound: Compound,
		supplier: str,
		catalogue: str,
		entry: str,
		amount: float,
		price: float,
		currency: str,
		purity: float,
		lead_time: int,
		smiles: str | None = None,
	) -> int | None:

		assert isinstance(compound, Compound), f'incompatible {compound=}'
		assert currency in ['GBP', 'EUR', 'USD'], f'incompatible {currency=}'

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
			self.execute(sql, (smiles, amount, supplier, catalogue, entry, lead_time, price, currency, purity, compound.id))

		except sqlite3.InterfaceError as e:
			logger.error(e)
			logger.debug((smiles, amount, supplier, catalogue, entry, lead_time, price, currency, purity, compound.id))
			raise

		except Exception as e:
			logger.exception(e)
			return None
			
		quote_id = self.cursor.lastrowid
		self.commit()
		return quote_id

	def insert_target(self,
		*,
		name: str
	) -> int:

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

		# except sqlite3.IntegrityError as e:
		# 	logger.warning(f"Target with {name=} already exists")

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

		payload = json.dumps(payload)

		self.update(table=table, id=id, key=f'{table}_metadata', value=payload, commit=commit)

	### SELECTION

	def select(self, query, table, multiple=False):

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

	def select_where(self, query, table, key, value, multiple=False, none='error', sort=None):

		if isinstance(value, str):
			value = f"'{value}'"

		if sort:
			sql = f'SELECT {query} FROM {table} WHERE {table}_{key}={value} ORDER BY {sort}'
		else:
			sql = f'SELECT {query} FROM {table} WHERE {table}_{key}={value}'

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
			logger.error(f'No entry in {table} with {table}_{key}={value}')
			return None

		return result

	def select_id_where(self, table, key, value, multiple=False):
		return self.select_where(query=f'{table}_id', table=table, key=key, value=value, multiple=multiple)

	def select_all_where(self, table, key, value, multiple=False):
		return self.select_where(query='*', table=table, key=key, value=value, multiple=multiple)

	### DELETION

	def delete_where(self, table, key, value):

		if isinstance(value, str):
			value = f"'{value}'"

		sql = f'DELETE FROM {table} WHERE {table}_{key}={value}'

		try:
			result = self.execute(sql)
			self.commit()

		except sqlite3.OperationalError as e:
			logger.var('sql',sql)
			raise

		return None


	### UPDATE

	def update(self, *, 
		table: str, 
		id: int, 
		key: str, 
		value, 
		commit: bool = True,
	):

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
		name: str | None = None,
		smiles: str | None = None,
	) -> Compound:
		
		if id is None:
			id = self.get_compound_id(name=name, smiles=smiles)

		if not id:
			logger.error(f'Invalid {id=}')
			return None

		query = 'compound_id, compound_name, compound_smiles, compound_base'
		entry = self.select_where(query=query, table=table, key='id', value=id)
		compound = Compound(self, *entry, metadata=None, mol=None)
		return compound

	def get_compound_id(self, 
		*,
		table: str = 'compound',
		name: str | None = None, 
		smiles: str | None = None, 
		similar: str | Compound | int | None = None, 
		threshold: float = 1.0,
	) -> int:

		if name:
			entry = self.select_id_where(table=table, key='name', value=name)

		if smiles:
			entry = self.select_id_where(table=table, key='smiles', value=smiles)

		if similar:
			raise NotImplementedError

		if entry:
			return entry[0]

		return None

	def get_pose(self,
		*,
		table: str = 'pose',
		id: int | None = None,
		name: str = None,
		longname: str | None = None,
		# smiles: str | None = None,
	) -> Pose:
		
		if id is None:
			id = self.get_pose_id(longname=longname, name=name)

		if not id:
			logger.error(f'Invalid {id=}')
			return None

		query = 'pose_id, pose_name, pose_longname, pose_smiles, pose_reference, pose_path, pose_compound, pose_target, pose_mol, pose_fingerprint'
		entry = self.select_where(query=query, table=table, key='id', value=id)
		pose = Pose(self, *entry)
		return pose

	def get_pose_id(self, 
		*,
		table: str = 'pose',
		longname: str | None = None, 
		name: str | None = None, 
	) -> int:

		if name:
			entry = self.select_id_where(table=table, key='name', value=name)

		elif longname:
			entry = self.select_id_where(table=table, key='longname', value=longname)
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
		
		query = 'quote_compound, quote_supplier, quote_catalogue, quote_entry, quote_amount, quote_price, quote_currency, quote_lead_time, quote_purity, quote_date, quote_smiles '
		entry = self.select_where(query=query, table=table, key='id', value=id, none=none)

		return Quote(
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
		return Target(db=self,id=id,name=self.get_target_name(id=id))

	def get_target_name(self,
		*,
		id:int,
	) -> str:

		table = 'target'
		payload, = self.select_where(query=f'{table}_name', table=table, key=f'id', value=id)
		return payload

	def get_target_id(self, 
		*,
		name: str, 
	) -> int:
			
		table = 'target'
		entry = self.select_id_where(table=table, key='name', value=name)

		if entry:
			return entry[0]

		return None

	def get_feature(self, 
		*,
		id=int
	) -> Feature:
		entry = self.select_all_where(table='feature', key='id', value=id)
		return Feature(*entry)

	### COMPOUND QUERY

	def query_substructure(self, query, fast=True, none='error'):

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

		from .cset import CompoundSubset

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
			cset = CompoundSubset(self, 'compound', ids)
			return cset, similarities

		ids = [r for r, in result]
		cset = CompoundSubset(self, 'compound', ids)

		return cset

	def query_exact(self, query):
		return self.query_similarity(query, 0.989, return_similarity=False)

	### COUNTING

	def count(self, table):

		sql = f"""SELECT COUNT(1) FROM {table}; """
		self.execute(sql)
		return self.cursor.fetchone()[0]

	def count_where(self, table, key, value):
		sql = f"""SELECT COUNT(1) FROM {table} WHERE {table}_{key} is {value}; """
		self.execute(sql)
		return self.cursor.fetchone()[0]

	### PRUNING

	def prune_reactions(self, compound, reactions):

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

	### PRINTING

	def list_tables(self):
		self.execute("SELECT name FROM sqlite_master WHERE type='table';")
		pprint(self.cursor.fetchall())

	def print_table(self, name):
		self.execute(f"SELECT * FROM {name}")
		pprint(self.cursor.fetchall())

	### DUNDERS
