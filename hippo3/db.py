
import mout
import mcol
import sqlite3
from sqlite3 import Error
from pprint import pprint
import json

from .compound import Compound
from .pose import Pose
from .reaction import Reaction

from pathlib import Path

import logging
logger = logging.getLogger('HIPPO')

class Database:

	def __init__(self, path: Path):

		assert isinstance(path, Path)
		
		logger.debug('hippo3.Database.__init__()')

		self._path = path
		self._connection = None
		self._cursor = None

		logger.debug(f'Database.path = {self.path}')


		try:
			# create connection
			path = path.resolve(strict=True)

		except FileNotFoundError:
			# doesn't exist
			# create a blank database
			self.connect()
			self.create_blank_db()

		else:
			# exists
			self.connect()
		
	### FACTORIES

	@classmethod
	def from_file(cls, path):

		logger.debug('hippo3.Database.from_file()')

		# load DB from file

		self = cls.__new__(cls)

		self.__init__()
		
		return self

	### PROPERTIES

	@property
	def path(self):
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

	def close(self):
		logger.debug('hippo3.Database.close()')
		if self.connection:
			self.connection.close()
		mout.success(f'Closed connection to @ {mcol.file}{self.path}')

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
		
		finally:
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
		self.connection.commit()

	### CREATION

	def create_blank_db(self):

		logger.out('Creating blank database...')
		self.create_table_compound()
		self.create_table_inspiration()
		self.create_table_reaction()
		self.create_table_reactant()
		self.create_table_pose()
		self.create_table_tag()
		self.create_table_quote()

	def create_table_compound(self):
		logger.debug('HIPPO.Database.create_table_compound()')

		sql = """CREATE TABLE compound(
			compound_id INTEGER PRIMARY KEY,
			compound_name TEXT,
			compound_smiles TEXT,
			compound_base INTEGER,
			compound_mol MOL,
			FOREIGN KEY (compound_base) REFERENCES compound(compound_id),
			CONSTRAINT UC_compound_name UNIQUE (compound_name),
			CONSTRAINT UC_compound_smiles UNIQUE (compound_smiles)
		);
		"""

		self.execute(sql)

		# if logging.root.level >= logging.DEBUG:
		# 	self.execute(f'PRAGMA table_info(compound);')
		# 	pprint(self.cursor.fetchall())

	def create_table_inspiration(self):
		logger.debug('HIPPO.Database.create_table_inspiration()')

		sql = """CREATE TABLE inspiration(
			inspiration_compound INTEGER,
			inspiration_pose INTEGER,
			FOREIGN KEY (inspiration_compound) REFERENCES compound(compound_id),
			FOREIGN KEY (inspiration_pose) REFERENCES pose(pose_id),
			CONSTRAINT UC_inspiration UNIQUE (inspiration_compound, inspiration_pose)
		);
		"""

		self.execute(sql)

		# if logging.root.level >= logging.DEBUG:
		# 	self.execute(f'PRAGMA table_info(inspiration);')
		# 	pprint(self.cursor.fetchall())

	def create_table_reaction(self):
		logger.debug('HIPPO.Database.create_table_reaction()')

		sql = """CREATE TABLE reaction(
			reaction_id INTEGER PRIMARY KEY,
			reaction_type TEXT,
			reaction_product INTEGER,
			reaction_product_amount REAL,
			FOREIGN KEY (reaction_product) REFERENCES compound(compound_id)
		);
		"""

		self.execute(sql)

		# if logging.root.level >= logging.DEBUG:
		# 	self.execute(f'PRAGMA table_info(reaction);')
		# 	pprint(self.cursor.fetchall())

	def create_table_reactant(self):
		logger.debug('HIPPO.Database.create_table_reactant()')

		sql = """CREATE TABLE reactant(
			reactant_amount REAL,
			reactant_reaction INTEGER,
			reactant_compound INTEGER,
			FOREIGN KEY (reactant_reaction) REFERENCES reaction(reaction_id)
			FOREIGN KEY (reactant_compound) REFERENCES compound(compound_id)
			CONSTRAINT UC_reactant UNIQUE (reactant_reaction, reactant_compound)
		);
		"""

		self.execute(sql)

		# if logging.root.level >= logging.DEBUG:
		# 	self.execute(f'PRAGMA table_info(reactant);')
		# 	pprint(self.cursor.fetchall())

	def create_table_pose(self):
		logger.debug('HIPPO.Database.create_table_pose()')

		sql = """CREATE TABLE pose(
			pose_id INTEGER PRIMARY KEY,
			pose_name TEXT,
			pose_longname TEXT,
			pose_smiles TEXT,
			pose_compound INTEGER,
			pose_path TEXT,
			pose_target TEXT,
			pose_mol BLOB,
			pose_fingerprint BLOB,
			FOREIGN KEY (pose_compound) REFERENCES compound(compound_id),
			CONSTRAINT UC_pose_longname UNIQUE (pose_longname)
		);
		"""

		self.execute(sql)

		# if logging.root.level >= logging.DEBUG:
		# 	self.execute(f'PRAGMA table_info(pose);')
		# 	pprint(self.cursor.fetchall())

	def create_table_tag(self):
		logger.debug('HIPPO.Database.create_table_tag()')

		sql = """CREATE TABLE tag(
			tag_name TEXT,
			tag_compound INTEGER,
			tag_pose INTEGER,
			FOREIGN KEY (tag_compound) REFERENCES compound(compound_id)
			FOREIGN KEY (tag_pose) REFERENCES pose(pose_id)
		);
		"""

		self.execute(sql)

		# if logging.root.level >= logging.DEBUG:
		# 	self.execute(f'PRAGMA table_info(tag);')
		# 	pprint(self.cursor.fetchall())

	def create_table_quote(self):
		logger.debug('HIPPO.Database.create_table_quote()')

		sql = """CREATE TABLE quote(
			quote_id INTEGER PRIMARY KEY,
			quote_amount REAL,
			quote_supplier TEXT,
			quote_catalogue TEXT,
			quote_entry TEXT,
			quote_lead_time INTEGER,
			quote_price REAL,
			quote_currency TEXT,
			quote_date TEXT,
			quote_compound INTEGER,
			FOREIGN KEY (quote_compound) REFERENCES compound(compound_id)
		);
		"""

		self.execute(sql)

		# if logging.root.level >= logging.DEBUG:
		# 	self.execute(f'PRAGMA table_info(quote);')
		# 	pprint(self.cursor.fetchall())

	### INSERTION

	def insert_compound(self, 
		name: str, 
		smiles: str, 
		base: Compound | None = None, 
		tags: None | list = None, 
		# inspirations=None,
	) -> int:

		# process the base
		assert isinstance(base, Compound) or base is None, f'incompatible base={base}'

		if base:
			base = base.id

		sql = """
		INSERT INTO compound(compound_name, compound_smiles, compound_base, compound_mol)
		VALUES(?1, ?2, ?3, mol_from_smiles(?2))
		"""

		try:
			self.execute(sql, (name, smiles, base))

		except sqlite3.IntegrityError as e:
			if 'UNIQUE constraint failed: compound.compound_name' in str(e):
				logger.error(f"Can't add compound with existing name \"{name}\"")
			elif 'UNIQUE constraint failed: compound.compound_smiles' in str(e):
				logger.error(f"Can't add compound with existing smiles \"{smiles}\"")
			else:
				logger.exception(e)

		except Exception as e:
			logger.exception(e)

		finally:
			
			compound_id = self.cursor.lastrowid

			if tags:
				for tag in tags:
					self.insert_tag(name=tag, compound=compound_id)

			return compound_id

		return None

	def insert_pose(self,
		compound: Compound,
		name: str,
		# smiles: str,
		target: str,
		path: str | None = None,
		tags: None | list = None,
	):

		assert isinstance(compound, Compound), f'incompatible {compound}'
		assert name, f'incompatible name={name}'
		# assert smiles

		longname = f'{target}_{compound.name}_{name}'

		if path:
			path = Path(path)
			path = path.resolve(strict=True)
			path = str(path)

		sql = """
		INSERT INTO pose(pose_name, pose_longname, pose_smiles, pose_compound, pose_target, pose_path)
		VALUES(?1, ?2, ?3, ?4, ?5, ?6)
		"""

		try:
			self.execute(sql, (name, longname, compound.smiles, compound.id, target, path))

		except sqlite3.IntegrityError as e:
			logger.error(f"Can't add pose with existing longname \"{longname}\"")

		except Exception as e:
			logger.exception(e)

		finally:
			
			pose_id = self.cursor.lastrowid

			if tags:
				for tag in tags:
					self.insert_tag(name=tag, pose=pose_id)

			return pose_id

		return None

	def insert_tag(self, 
		name, 
		compound: int = None, 
		pose: int = None
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
				logger.error(f"Can't add compound with existing name \"{name}\"")
			elif 'UNIQUE constraint failed: compound.compound_smiles' in str(e):
				logger.error(f"Can't add compound with existing smiles \"{smiles}\"")
			else:
				logger.exception(e)

		except Exception as e:
			logger.exception(e)		

		finally:
			tag_id = self.cursor.lastrowid
			return tag_id

		return None

	def insert_inspiration(self, 
		compound: Compound, 
		pose: Pose,
	) -> int:

		assert isinstance(compound, Compound), f'incompatible compound={compound}'
		assert isinstance(pose, Pose), f'incompatible pose={pose}'

		# assert isinstance(compound_id, int) or compound_id is None
		# assert isinstance(pose_id, int) or pose_id is None

		sql = """
		INSERT INTO inspiration(inspiration_compound, inspiration_pose)
		VALUES(?1, ?2)
		"""

		try:
			self.execute(sql, (compound.id, pose.id))

		except sqlite3.IntegrityError as e:
			logger.warning(f"Can't add existing inspiration: {compound} {pose}")

		except Exception as e:
			logger.exception(e)

		finally:
			inspiration_id = self.cursor.lastrowid
			return inspiration_id

		return None

	def insert_reaction(self,
		type: str,
		product: Compound,
		product_amount: float,
	) -> int:

		assert isinstance(product, Compound), f'incompatible {product=}'
		assert isinstance(type, str), f'incompatible {type=}'

		sql = """
		INSERT INTO reaction(reaction_type, reaction_product, reaction_product_amount)
		VALUES(?1, ?2, ?3)
		"""

		try:
			self.execute(sql, (type, product.id, product_amount))

		except Exception as e:
			logger.exception(e)

		finally:
			
			reaction_id = self.cursor.lastrowid
			return reaction_id

		return None

	def insert_reactant(self,
		compound: Compound,
		reaction: Reaction,
		amount: float,
	) -> int:

		assert isinstance(compound, Compound), f'incompatible {compound=}'
		assert isinstance(reaction, Reaction), f'incompatible {reaction=}'

		sql = """
		INSERT INTO reactant(reactant_amount, reactant_reaction, reactant_compound)
		VALUES(?1, ?2, ?3)
		"""

		try:
			self.execute(sql, (amount, reaction.id, compound.id))

		except Exception as e:
			logger.exception(e)

		finally:
			
			reactant_id = self.cursor.lastrowid
			return reactant_id

		return None

	def insert_quote(self,
		compound: Compound,
		supplier: str,
		catalogue: str,
		entry: str,
		amount: float,
		price: float,
		currency: str,
		lead_time: int,
	):

		assert isinstance(compound, Compound), f'incompatible {compound=}'
		assert currency in ['GBP', 'EUR', 'USD'], f'incompatible {currency=}'
		
		sql = """
		INSERT INTO quote(
			quote_amount,
			quote_supplier,
			quote_catalogue,
			quote_entry,
			quote_lead_time,
			quote_price,
			quote_currency,
			quote_compound,
			quote_date
		)
		VALUES(?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, date());
		"""

		try:
			self.execute(sql, (amount, supplier, catalogue, entry, lead_time, price, currency, compound.id))

		except Exception as e:
			logger.exception(e)

		finally:
			
			quote_id = self.cursor.lastrowid
			return quote_id

		return None

	### SELECTION

	def select_where(self, query, table, key, value, multiple=False, none='error'):
		if isinstance(value, str):
			value = f"'{value}'"

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

	### GETTERS

	def get_compound(self,
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

		query = 'compound_id, compound_name, compound_smiles, compound_base, mol_to_binary_mol(compound_mol)'
		entry = self.select_where(query=query, table=table, key='id', value=id)
		compound = Compound(self, *entry)
		return compound

	def get_compound_id(self, 
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
		table: str = 'pose',
		id: int | None = None,
		longname: str | None = None,
		smiles: str | None = None,
	) -> Pose:
		
		if id is None:
			id = self.get_pose_id(longname=longname)

		if not id:
			logger.error(f'Invalid {id=}')
			return None

		query = 'pose_id, pose_name, pose_longname, pose_smiles, pose_compound, pose_target, pose_path, pose_mol, pose_fingerprint'
		entry = self.select_where(query=query, table=table, key='id', value=id)
		pose = Pose(self, *entry)
		return pose

	def get_pose_id(self, 
		table: str = 'pose',
		longname: str | None = None, 
	) -> int:
			
		entry = self.select_id_where(table=table, key='longname', value=longname)

		if entry:
			return entry[0]

		return None

	def get_reaction(self,
		table: str = 'reaction',
		id: int | None = None,
	) -> Reaction:
		
		if not id:
			logger.error(f'Invalid {id=}')
			return None

		query = 'reaction_id, reaction_type, reaction_product, reaction_product_amount'
		entry = self.select_where(query=query, table=table, key='id', value=id)
		reaction = Reaction(self, *entry)
		return reaction

	def get_quote(self,
		table: str = 'quote',
		id: int | None = None,
	) -> list[dict]:
		
		query = 'quote_compound, quote_supplier, quote_catalogue, quote_entry, quote_amount, quote_price, quote_currency, quote_lead_time, quote_date'
		entry = self.select_where(query=query, table=table, key='id', value=id)

		return dict(
			compound=entry[0],
			supplier=entry[1],
			catalogue=entry[2],
			entry=entry[3],
			amount=entry[4],
			price=entry[5],
			currency=entry[6],
			lead_time=entry[7],
			date=entry[8],
		)

		return entry

	### COUNTING

	def count(self, table):

		sql = f"""SELECT COUNT(1) FROM {table}; """
		self.execute(sql)
		return self.cursor.fetchone()[0]

	### PRINTING

	def list_tables(self):
		self.execute("SELECT name FROM sqlite_master WHERE type='table';")
		pprint(self.cursor.fetchall())

	def print_table(self, name):
		self.execute(f"SELECT * FROM {name}")
		pprint(self.cursor.fetchall())

	### DUNDERS
