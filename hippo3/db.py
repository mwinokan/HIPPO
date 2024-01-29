
import mout
import mcol
import sqlite3
from sqlite3 import Error
from pprint import pprint

class Database:

	def __init__(self, path):
		mout.debug('Database.__init__()')

		self._path = path
		self._connection = None
		self._cursor = None

		mout.var('Database.path', self.path, valCol=mcol.file)

		# create connection
		self._connect()

		# create a blank database
		self._create_blank_db()
		
	### FACTORIES

	@classmethod
	def from_file(cls, path):

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
			self._connect()
		return self._connection
	
	@property
	def cursor(self):
		if not self._cursor:
			serlf._cursor = self.connection.cursor()
		return self._cursor

	### PUBLIC METHODS / API CALLS

	def close(self):
		if self.connection:
			self.connection.close()
		mout.success(f'Closed connection to @ {mcol.file}{self.path}')

	def register_compound(self,
			name,
			smiles,
			stereo_smiles,
			orig_smiles,
			mol=None,
		):

		# for field in TABLE_COMPOUND['fields']:
		# 	assert field in kwargs, field

		self._create_entry('COMPOUND', 
			name=name, 
			smiles=smiles, 
			stereo_smiles=stereo_smiles, 
			orig_smiles=orig_smiles,
			mol=mol,
		)

	### INTERNAL METHODS

	def _connect(self):
		mout.debug('Database._connect()')

		conn = None

		try:
			conn = sqlite3.connect(self.path)
			mout.debug(f'{sqlite3.version=}')

			conn.enable_load_extension(True)
			# conn.load_extension('chemicalite')
			conn.load_extension('/Users/tfb64483/Software/miniconda3/envs/py310/lib/chemicalite.dylib')
			conn.enable_load_extension(False)

			mout.success(f'Database connected @ {mcol.file}{self.path}')
		
		except Error as e:
			mout.error(e)
		
		finally:
			self._connection = conn
			self._cursor = conn.cursor()

	def _execute(self, sql, payload=None, debug=True):
		try:
			if payload:
				if debug:
					mout.debug(f'SQL: {sql}')
					mout.debug(f'payload: {payload}')
				return self.cursor.execute(sql, payload)
			else:
				if debug:
					mout.debug(f'SQL: {sql}')
				return self.cursor.execute(sql)
		except Error as e:
			mout.error(str(e))

	def _commit(self):
		self.connection.commit()

	def _list_tables(self):
		self._execute("SELECT name FROM sqlite_master WHERE type='table';")
		pprint(self._cursor.fetchall())

	def _print_table(self, name):
		self._execute(f"SELECT * FROM {name}")
		pprint(self._cursor.fetchall())

	def _create_blank_db(self):

		# create empty tables
		self._create_table(TABLE_COMPOUND)

	def _create_table(self, template):

		sql = []

		n = len(template['fields'])

		sql.append(f'CREATE TABLE IF NOT EXISTS {template["name"]} (')
		sql.append(f'id integer PRIMARY KEY,')

		for i, (key, value) in enumerate(template['fields'].items()):

			mout.debug(f'{i}, {key}, {value}')

			line = []

			line.append(key)
			line.append(value)

			if i != n-1:
				line.append(',')
			
			line = ' '.join(line)

			sql.append(line)
		
		sql.append(');')

		sql = '\n'.join(sql)

		mout.debug(sql)

		self._execute(sql)

		self._execute(f'PRAGMA table_info({template["name"]});')
		pprint(self._cursor.fetchall())

	def _create_entry(self, table, **kwargs):

		sql = []

		n = len(kwargs)

		sql.append(f'INSERT INTO {table}(')

		for i, key in enumerate(kwargs):

			mout.debug(f'{i}, {key}')

			line = []

			line.append(key)

			if i != n-1:
				line.append(',')
			
			line = ' '.join(line)

			sql.append(line)
		
		sql.append(')')

		sql.append('VALUES(')

		sql.append(','.join(['?' for _ in range(n)]))

		sql.append(')')

		sql = '\n'.join(sql)

		mout.debug(sql)
		mout.debug(list(kwargs.values()))

		self._execute(sql, list(kwargs.values()))

		return self.cursor.lastrowid

	### DUNDERS

TABLE_COMPOUND = {
	'name': 'COMPOUND',
	'fields': {
		'name': 'TEXT NOT NULL',
		'smiles': 'TEXT NOT NULL',
		'stereo_smiles': 'TEXT NOT NULL',
		'orig_smiles': 'TEXT NOT NULL',
		'mol': 'MOL', # Use Chemicalite
		# 'INT': 'fp_1024', # Use Chemicalite
		# 'INT': 'canonical_smiles_hash',
		# 'REF': {'COMPOUND': 'parent'},
		# 'REFS': {'TAG': 'tags'},
		# 'REFS': {'POSE': 'poses'},
		# 'REFS': {'POSE': 'inspirations'},
		# 'REFS': {'PURCHASE': 'purchase_data'},
		# 'REFS': {'REACTION': 'reactions'},
	}
}
