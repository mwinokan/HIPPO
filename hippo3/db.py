
import mout
import mcol
import sqlite3
from sqlite3 import Error
from pprint import pprint
import json
from .compound import Compound

import logging
logger = logging.getLogger('HIPPO')

class Database:

	def __init__(self, path):
		
		logger.debug('hippo3.Database.__init__()')

		self._path = path
		self._connection = None
		self._cursor = None

		logger.debug(f'Database.path = {self.path}')

		# create connection
		self.connect()

		# create a blank database
		self.create_blank_db()
		
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

		self.create_table_compound()
		self.create_table_inspiration()
		self.create_table_reaction()
		self.create_table_reactant()
		self.create_table_pose()
		self.create_table_tag()
		self.create_table_purchase()

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

		if logging.root.level >= logging.DEBUG:
			self.execute(f'PRAGMA table_info(compound);')
			pprint(self.cursor.fetchall())

	def create_table_inspiration(self):
		logger.debug('HIPPO.Database.create_table_inspiration()')

		sql = """CREATE TABLE inspiration(
			inspiration_compound INTEGER,
			inspiration_pose INTEGER,
			FOREIGN KEY (inspiration_compound) REFERENCES compound(compound_id)
			FOREIGN KEY (inspiration_pose) REFERENCES pose(pose_id)
		);
		"""

		self.execute(sql)

		if logging.root.level >= logging.DEBUG:
			self.execute(f'PRAGMA table_info(inspiration);')
			pprint(self.cursor.fetchall())

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

		if logging.root.level >= logging.DEBUG:
			self.execute(f'PRAGMA table_info(reaction);')
			pprint(self.cursor.fetchall())

	def create_table_reactant(self):
		logger.debug('HIPPO.Database.create_table_reactant()')

		sql = """CREATE TABLE reactant(
			reactant_amount REAL,
			reactant_reaction INTEGER,
			reactant_compound INTEGER,
			FOREIGN KEY (reactant_reaction) REFERENCES reaction(reaction_id)
			FOREIGN KEY (reactant_compound) REFERENCES compound(compound_id)
		);
		"""

		self.execute(sql)

		if logging.root.level >= logging.DEBUG:
			self.execute(f'PRAGMA table_info(reactant);')
			pprint(self.cursor.fetchall())

	def create_table_pose(self):
		logger.debug('HIPPO.Database.create_table_pose()')

		sql = """CREATE TABLE pose(
			pose_id INTEGER PRIMARY KEY,
			pose_name TEXT,
			pose_smiles TEXT,
			pose_compound INTEGER,
			pose_path TEXT,
			pose_mol BLOB,
			pose_fingerprint BLOB,
			FOREIGN KEY (pose_compound) REFERENCES compound(compound_id)
		);
		"""

		self.execute(sql)

		if logging.root.level >= logging.DEBUG:
			self.execute(f'PRAGMA table_info(pose);')
			pprint(self.cursor.fetchall())

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

		if logging.root.level >= logging.DEBUG:
			self.execute(f'PRAGMA table_info(tag);')
			pprint(self.cursor.fetchall())

	def create_table_purchase(self):
		logger.debug('HIPPO.Database.create_table_purchase()')

		sql = """CREATE TABLE purchase(
			purchase_amount REAL,
			purchase_supplier TEXT,
			purchase_catalogue TEXT,
			purchase_entry TEXT,
			purchase_lead_time INTEGER,
			purchase_price REAL,

			purchase_compound INTEGER,
			FOREIGN KEY (purchase_compound) REFERENCES compound(compound_id)
		);
		"""

		self.execute(sql)

		if logging.root.level >= logging.DEBUG:
			self.execute(f'PRAGMA table_info(purchase);')
			pprint(self.cursor.fetchall())

	### INSERTION

	def add_compound(self, 
		name, 
		smiles, 
		base=None, 
		tags=None, 
		# inspirations=None,
	):

		# process the base
		assert isinstance(base, int) or base is None

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
					self.add_tag(name=tag, compound=compound_id)

			# if inspirations:
			# 	for inspiration in inspirations:
			# 		self.add_inspiration()

			return compound_id

		return None

	def add_pose(self,
		name: str,
		smiles: str,
		compound: int,
		path: str | None = None,
	):

		assert name
		assert smiles
		assert compound

		sql = """
		INSERT INTO pose(pose_name, pose_smiles, pose_compound, pose_path)
		VALUES(?1, ?2, ?3, ?4)
		"""

		try:
			self.execute(sql, (name, smiles, compound, path))

		except Exception as e:
			logger.exception(e)

		finally:
			
			pose_id = self.cursor.lastrowid
			return pose_id

		return None

	def add_tag(self, name, compound=None, pose=None):	
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

	def add_inspiration(self, compound_id, pose_id):

		assert isinstance(compound_id, int) or compound_id is None
		assert isinstance(pose_id, int) or pose_id is None

		sql = """
		INSERT INTO inspiration(inspiration_compound, inspiration_pose)
		VALUES(?1, ?2)
		"""

		try:
			self.execute(sql, (compound_id, pose_id))

		except Exception as e:
			logger.exception(e)

		finally:
			inspiration_id = self.cursor.lastrowid
			return inspiration_id

		return None

	### SELECTION

	def select_where(self, query, table, key, value, multiple=False):
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

		if not result:
			logger.error(f'No entry in {table} with {table}_{key}={value}')
			return None

		return result

	def select_id_where(self, table, key, value, multiple=False):
		return self.select_where(query=f'{table}_id', table=table, key=key, value=value, multiple=multiple)

	def select_all_where(self, table, key, value, multiple=False):
		return self.select_where(query='*', table=table, key=key, value=value, multiple=multiple)

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
