
from .tools import df_row_to_dict

from .compound import Compound
from .db import Database

import os

import logging
logger = logging.getLogger('HIPPO')

class CompoundSet:

	# class to access entries in database tables containing compounds

	def __init__(self, 
		db: Database, 
		table: str = 'compound',
	) -> None:
		
		self._db = db
		self._table = table
		
	### FACTORIES

	### PROPERTIES

	@property
	def db(self) -> Database:
		return self._db
	
	@property
	def table(self) -> str:
		return self._table

	### METHODS

	def get_compound(self,
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
		entry = self.db.select_where(query=query, table='compound', key='id', value=id)
		compound = Compound.from_db_entry(self, *entry)
		return compound

	def get_compound_id(self, 
		name: str | None = None, 
		smiles: str | None = None, 
		similar: str | Compound | int | None = None, 
		threshold: float = 1.0,
	) -> int:

		if name:
			entry = self.db.select_id_where(table=self.table, key='name', value=name)

		if smiles:
			entry = self.db.select_id_where(table=self.table, key='smiles', value=smiles)

		if similar:
			raise NotImplementedError

		if entry:
			return entry[0]

		return None

	### DUNDERS

	def __getitem__(self, key) -> Compound:
		
		match key:

			case int():
				return self.get_compound(id=key)

			case str():
				return self.get_compound(name=key)

			case _:
				logger.error(f'Unsupported type for CompoundSet.__getitem__(): {type(key)}')

		return None

	def __repr__(self) -> str:
		return f'CompoundSet("compound")'

	def __len__(self) -> int:
		return self.db.count(self.table)
