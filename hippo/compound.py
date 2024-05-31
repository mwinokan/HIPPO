
import mcol

from rdkit import Chem

from .pose import Pose
# from .pset import PoseSet
from .tags import TagSet
# from .rset import ReactionSet
from .target import Target
from .quote import Quote

import logging
logger = logging.getLogger('HIPPO')


class Compound:

	"""A :class:`.Compound` represents a ligand/small molecule with stereochemistry removed and no atomic coordinates. I.e. it represents the chemical structure. It's name is always an InChiKey. If a compound is an elaboration it can have a :meth:`.Compound.base` property which is another :class:`.Compound`. :class:`.Compound` objects are target-agnostic and can be linked to any number of catalogue entries (:class:`.Quote`) or synthetic pathways (:class:`.Reaction`).

	:class:`.Compound` objects should not be created directly. Instead use :meth:`.HIPPO.register_compound` or :meth:`.HIPPO.compounds`
	"""

	def __init__(self,
			animal,
			db,
			id: int,
			inchikey: str,
			alias: str,
			smiles: str,
			base: int,
			mol: Chem.Mol | bytes | None = None,
			metadata: dict | None = None,
	):
		
		# from compound table
		self._id = id
		self._inchikey = inchikey
		self._alias = alias
		self._smiles = smiles
		self._animal = animal
		self._base_id = base
		self._base = None
		self._alias = alias
		self._tags = None
		self._table = 'compound'
		self._metadata = metadata
		
		# computed properties
		self._num_heavy_atoms = None
		self._num_rings = None
		self._formula = None
		self._molecular_weight = None
		
		if isinstance(mol, bytes):
			mol = Chem.Mol(mol)
			
		self._mol = mol

		self._db = db
		
	### FACTORIES

	### PROPERTIES

	@property
	def id(self) -> int:
		"""Returns the compound's database ID"""
		return self._id
	
	@property
	def inchikey(self) -> str:
		"""Returns the compound's InChiKey"""
		return self._inchikey
	
	@property
	def name(self) -> str:
		"""Returns the compound's InChiKey"""
		return self._inchikey
	
	@property
	def smiles(self) -> str:
		"""Returns the compound's (flattened) smiles"""
		return self._smiles

	@property
	def alias(self) -> str:
		"""Returns the compound's alias"""
		return self._alias

	@alias.setter
	def alias(self, alias: str) -> None:
		"""Set the compound's alias"""
		assert isinstance(alias, str)
		self._alias = alias
		self.db.update(table='compound', id=self.id, key='compound_alias', value=alias, commit=commit)
	
	@property
	def mol(self) -> Chem.Mol:
		"""Returns the compound's RDKit Molecule"""
		if self._mol is None:
			mol, = self.db.select_where(query='mol_to_binary_mol(compound_mol)', table='compound', key='id', value=self.id, multiple=False)
			self._mol = Chem.Mol(mol)
		return self._mol

	@property
	def num_heavy_atoms(self):
		"""Get the number of heavy atoms"""
		if self._num_heavy_atoms is None:
			self._num_heavy_atoms = self.db.get_compound_computed_property('num_heavy_atoms', self.id)
		return self._num_heavy_atoms

	@property
	def molecular_weight(self):
		"""Get the molecular weight"""
		if self._molecular_weight is None:
			self._molecular_weight = self.db.get_compound_computed_property('molecular_weight', self.id)
		return self._molecular_weight

	@property
	def num_rings(self):
		"""Get the number of rings"""
		if self._num_rings is None:
			self._num_rings = self.db.get_compound_computed_property('num_rings', self.id)
		return self._num_rings

	@property
	def formula(self):
		"""Get the chemical formula"""
		if self._formula is None:
			self._formula = self.db.get_compound_computed_property('formula', self.id)
		return self._formula

	@property
	def atomtype_dict(self):
		from .tools import formula_to_atomtype_dict
		return formula_to_atomtype_dict(self.formula)

	@property
	def num_atoms_added(self):
		"""Calculate the number of atoms added relative to the base compound"""
		assert (b_id := self._base_id), f'{self} has no base defined'
		n_e = self.num_heavy_atoms
		n_b = self.db.get_compound_computed_property('num_heavy_atoms', b_id)
		return n_e - n_b
	
	@property
	def metadata(self) -> dict:
		"""Returns the compound's metadata dict"""
		if self._metadata is None:
			self._metadata = self.db.get_metadata(table='compound', id=self.id)
		return self._metadata

	@property
	def db(self):
		"""Returns a pointer to the parent database"""
		return self._db

	@property
	def tags(self) -> TagSet:
		"""Returns the compound's tags"""
		if not self._tags:
			self._tags = self.get_tags()
		return self._tags

	@property
	def poses(self):
		"""Returns the compound's poses"""
		return self.get_poses()

	@property
	def best_placed_pose(self):
		return self.poses.best_placed_pose

	@property
	def num_poses(self) -> int:
		"""Returns the number of associated poses"""
		return self.db.count_where(table='pose', key='compound', value=self.id)

	@property
	def num_reactions(self) -> int:
		"""Returns the number of associated reactions (product)"""
		return self.db.count_where(table='reaction', key='product', value=self.id)

	@property
	def num_reactant(self) -> int:
		"""Returns the number of associated reactions (reactant)"""
		return self.db.count_where(table='reactant', key='compound', value=self.id)

	@property
	def base(self):
		"""Returns the base compound for this elaboration"""
		if self._base_id and self._base is None:
			self._base = self.db.get_compound(id=self._base_id)
		return self._base

	@base.setter
	def base(self, b):
		"""Set the base compound for this elaboration"""
		self.set_base(b)

	@property
	def reactions(self):
		"""Returns the reactions resulting in this compound"""
		return self.get_reactions(none=False)

	@property
	def dict(self) -> dict:
		"""Returns a dictionary of this compound"""
		return self.get_dict()

	@property
	def table(self):
		return self._table

	@property
	def elabs(self):
		"""Return the set of elaborations based on this compound"""
		ids = self.db.select_where(query='compound_id', table='compound', key='base', value=self.id, multiple=True)
		ids = [q for q, in ids]
		from .cset import CompoundSet
		return CompoundSet(self.db, ids)

	@property
	def is_base(self):
		"""Is this Compound the basis for any elaborations?"""
		return bool(self.db.select_where(query='compound_id', table='compound', key='base', value=self.id, multiple=False, none='quiet'))

	@property
	def is_elab(self):
		"""Is this Compound the based on any other compound?"""
		return bool(self.base)

	@property
	def is_product(self):
		"""Is this Compound a product of at least one reaction"""
		return bool(self.get_reactions(none=False))
	
	
	### METHODS

	# def get_dict(self, mol=False):
	# 	"""Returns a dictionary of this compound"""

	# 	serialisable_fields = ['id','inchikey','alias','smiles']

	# 	if mol:
	# 		serialisable_fields.append('mol')

	# 	data = {}
	# 	for key in serialisable_fields:
	# 		data[key] = getattr(self, key)

	# 	if base := self.base:
	# 		data['base'] = base.inchikey

	# 	if metadata := self.metadata:
	# 		for key in metadata:
	# 			data[key] = metadata[key]

	# 	return data

	def add_stock(self, 
		amount: float,
		*,
		purity: float | None = None,
		entry: str | None = None,
		location: str | None = None,
		return_quote: bool = True,
	) -> int | Quote:

		"""Register a certain quantity of compound stock in the Database."""

		assert amount

		# search for existing in stock quotes
		existing = self.get_quotes(supplier='Stock', df=False)

		# supersede old in stock records
		if existing:
			delete = set()
			not_deleted = 0
			for quote in existing:

				if any([
					quote.entry != entry,
					quote.purity != purity,
					quote.catalogue != location,
				]):
					not_deleted += 1
					continue

				delete.add(quote.id)

			delete_str = str(tuple(delete)).replace(',)',')')

			self.db.delete_where(table='quote', key=f'quote_id IN {delete_str}')

			if delete:
				logger.warning(f'Removed {len(delete)} existing In-Stock Quotes')
			
			if not_deleted:
				logger.warning(f'Did not remove {not_deleted} existing In-Stock Quotes with differing entry/purity/location')

		# insert the new quote
		quote_id = self.db.insert_quote(
			compound=self.id,
			price=0,
			lead_time=0,
			currency=None,
			supplier='Stock',
			catalogue=location,
			entry=entry,
			amount=amount,
			purity=purity,
		)

		if return_quote:
			self.db.get_quote(id=quote_id)
		else:
			return quote_id

	def get_tags(self) -> set:
		tags = self.db.select_where(query='tag_name', table='tag', key='compound', value=self.id, multiple=True, none='quiet')
		return TagSet(self, {t[0] for t in tags}, commit=False)

	def get_quotes(self, min_amount=None, supplier=None, max_lead_time=None, none='quiet', pick_cheapest=False, df=False) -> list[dict]:
		"""Get all quotes associated to this compound"""

		quote_ids = self.db.select_where(query='quote_id', table='quote', key='compound', value=self.id, multiple=True, none=none)

		if quote_ids:
			quotes = [self.db.get_quote(id=q[0]) for q in quote_ids]
		else:
			return None

		if supplier:
			quotes = [q for q in quotes if q.supplier == supplier]

		if max_lead_time:
			quotes = [q for q in quotes if q.lead_time <= max_lead_time]

		if min_amount:
			suitable_quotes = [q for q in quotes if q.amount >= min_amount]
			
			if not suitable_quotes:
				# logger.debug(f'No quote available with amount >= {min_amount} mg')
				quotes = [Quote.combination(min_amount, quotes)]

			else:
				quotes = suitable_quotes

		if pick_cheapest:
			return sorted(quotes, key=lambda x: x.price)[0]

		if df:
			from pandas import DataFrame
			return DataFrame([q.dict for q in quotes]).drop(columns='compound')
		
		return quotes

	def get_reactions(self, none='error', as_reactant=False, permitted_reactions=None) -> list:
		"""Get the associated reactions as product, unless as_reactant is True."""

		from .rset import ReactionSet

		if as_reactant:
			reaction_ids = self.db.select_where(query='reactant_reaction', table='reactant', key='compound', value=self.id, multiple=True, none=none)
		else:
			reaction_ids = self.db.select_where(query='reaction_id', table='reaction', key='product', value=self.id, multiple=True, none=none)

		reaction_ids = [q for q, in reaction_ids]

		if permitted_reactions:
			reaction_ids = [i for i in reaction_ids if i in permitted_reactions]

		return ReactionSet(self.db, reaction_ids)

	def get_poses(self, 
		target: str = None
	) -> list[Pose]:

		pose_ids = self.db.select_where(query='pose_id', table='pose', key='compound', value=self.id, multiple=True, none=False)

		# if not pose_ids:
			# return None

		# poses = [self.db.get_pose(id=q[0]) for q in pose_ids]

		from .pset import PoseSet

		return PoseSet(self.db, [q[0] for q in pose_ids])

	def get_dict(self, mol=True, reactions=False, metadata=True, count_by_target=False, poses=True):
		
		"""Returns a dictionary representing this Compound"""

		# serialisable_fields = ['id','alias', 'inchikey', 'smiles', 'num_poses', 'num_reactant', 'num_reactions']
		serialisable_fields = ['id','alias', 'inchikey', 'smiles', 'num_reactant', 'num_reactions']

		# poses
		# reactions
		# poses.targets

		# assert not reactions

		data = {}
		for key in serialisable_fields:
			data[key] = getattr(self, key)

		if mol:
			try:
				data['mol'] = self.mol
			except InvalidMolError:
				data['mol'] = None

		if self.base:
			data['base'] = self.base.name
		else:
			data['base'] = None

		data['tags'] = self.tags
		
		if poses:

			poses = self.poses
			
			if poses:

				# data['poses'] = [a if a else i for a,i in zip(poses.aliases, poses.inchikeys)]
				# data['poses'] = self.poses.names
				data['poses'] = poses.ids
				data['targets'] = poses.target_names
			
				if count_by_target:
					target_ids = poses.target_ids

					for target in self._animal.targets:
						t_poses = poses(target=target.id) or []
						data[f'#poses {target.name}'] = len(t_poses)
		
		if metadata and (metadict := self.metadata):
			for key in metadict:
				data[key] = metadict[key]

		return data

	def set_base(self, base, commit=True):
		if not isinstance(base, int):
			assert base._table == 'compound'
			base = base.id
		self._base_id = base
		self.db.update(table='compound', id=self.id, key='compound_base', value=base, commit=commit)

	def as_ingredient(self, amount, max_lead_time=None, supplier=None):
		"""Convert this compound into an Ingredient object"""
		
		quote = self.get_quotes(
			pick_cheapest=True, 
			min_amount=amount, 
			max_lead_time=max_lead_time, 
			supplier=supplier,
			none='quiet',
		)

		if not quote:
			quote = None

		# else:
			# quote = quote.id
		
		return Ingredient(self.db, self.id, amount, quote, supplier, max_lead_time)

	def draw(self, align_substructure: bool = False):
		"""Display this compound (and its base if it has one)"""

		if (base := self.base):

			# print(self.base)

			from molparse.rdkit import draw_mcs

			data = {self.base.smiles:f'{base} (base)', self.smiles:str(self)}

			if len(data) == 2:

				# print(data)

				drawing = draw_mcs(data, align_substructure=align_substructure, show_mcs=False, highlight_common=False)
				display(drawing)

			else:
				logger.error(f'Problem drawing {base.id=} vs {self.id=}, self referential?')
				display(self.mol)
		
		else:
			display(self.mol)

	def classify(self, draw=True):
		"""Find RDKit Fragments within the compound molecule and draw them (or just return a list of (descriptor, count) tuples)"""
		from molparse.rdkit import classify_mol
		return classify_mol(self.mol, draw=draw)

	def murcko_scaffold(self, generic=False):
		"""Get the rdkit MurckoScaffold for this compound"""

		from rdkit.Chem.Scaffolds import MurckoScaffold

		scaffold = MurckoScaffold.GetScaffoldForMol(self.mol)

		if generic:
			scaffold = MurckoScaffold.MakeScaffoldGeneric(scaffold)

		return scaffold

	def summary(self, metadata: bool = True, draw: bool = True):
		"""Print a summary of this compound"""
		logger.header(repr(self))

		logger.var('inchikey', self.inchikey)
		logger.var('alias', self.alias)
		logger.var('smiles', self.smiles)
		logger.var('base', self.base)
		
		logger.var('is_base', self.is_base)
		logger.var('num_heavy_atoms', self.num_heavy_atoms)
		logger.var('num_rings', self.num_rings)
		logger.var('formula', self.formula)

		poses = self.poses
		logger.var('#poses', len(poses))
		if poses:
			logger.var('targets', poses.targets)
		
		logger.var('#reactions (product)', self.num_reactions)
		logger.var('#reactions (reactant)', self.num_reactant)
		
		logger.var('tags', self.tags)
		
		if metadata:
			logger.var('metadata', str(self.metadata))
		
		if draw:
			return self.draw()

	def place(self,
		*,
		target: str | int | Target | None = None,
		inspirations: list[Pose] | None = None,
		reference: Pose | None = None,
		max_ddG: float = 0.0,
		max_RMSD: float = 2.0,
		output_dir: str = 'wictor_place',
		tags = None,
		metadata = None,
		overwrite = False,
	) -> Pose:
		"""Generate a new pose for this compound using Fragmenstein."""
		
		from fragmenstein import Monster, Wictor
		from pathlib import Path

		tags = tags or []
		metadata = metadata or {}

		# get required data
		smiles = self.smiles

		inspirations = inspirations or self.poses[0].inspirations
		target = target or self.poses[0].target.name
		reference = reference or self.poses[0].reference

		inspiration_mols = [c.mol for c in inspirations]
		protein_pdb_block = reference.protein_system.pdb_block_with_alt_sites
				
		# create the victor
		victor = Wictor(hits=inspiration_mols, pdb_block=protein_pdb_block)
		victor.work_path = output_dir
		victor.enable_stdout(logging.CRITICAL)

		# do the placement
		victor.place(smiles, long_name=self.name)

		# metadata
		metadata['ddG'] = victor.energy_score['bound']['total_score'] - victor.energy_score['unbound']['total_score']
		metadata['RMSD'] = victor.mrmsd.mrmsd

		# print(victor.energy_score)

		if metadata['ddG'] > max_ddG:
			return None

		if metadata['RMSD'] > max_RMSD:
			return None

		# register the pose
		pose = self._animal.register_pose(
			compound=self,
			target=target,
			path=Path(victor.work_path) / self.name / f'{self.name}.minimised.mol',
			inspirations = inspirations,
			reference=reference,
			tags=tags,
			metadata=metadata,	
		)

		if overwrite:
			ids = [p.id for p in self.poses if p.id != pose.id]
			for i in ids:
				self.db.delete_where(table='pose', key="id", value=i)
			logger.success(f'Successfully posed {self} (and deleted old poses)')
		else:
			logger.success(f'Successfully posed {self}')

		return pose

	### DUNDERS

	def __str__(self):
		return f'C{self.id}'

	def __repr__(self):
		return f'{mcol.bold}{mcol.underline}{self} "{self.name}"{mcol.unbold}{mcol.ununderline}'

	def __eq__(self, other):
		return self.id == other.id

class Ingredient:

	"""An ingredient is a :class:`.Compound` with a fixed quanitity and an attached quote.

	Create one from a :meth:`.Compound.as_ingredient`"""

	_table = 'ingredient'

	def __init__(self, db, compound, amount, quote, max_lead_time=None, supplier=None):

		assert compound

		self._db = db
			
		# don't store inherited compound in memory until needed
		self._compound = None
	
		if isinstance(compound, Compound):
			self._compound_id = compound.id
			self._compound = None
		else:
			self._compound_id = compound

		if isinstance(quote, Quote):

			if id := quote.id:
				self._quote_id = quote.id
				self._quote = None

			else:
				# logger.debug('Initialising Ingredient with estimated quote')
				self._quote_id = None
				self._quote = quote

		elif quote is None:
			self._quote_id = None
			self._quote = None
		
		else:
			self._quote_id = int(quote)
			self._quote = None
		
		# self._id = inherit.id
		# self._inchikey = inherit.inchikey
		# self._alias = inherit.alias
		# self._smiles = inherit.smiles
		# self._base = inherit.base			
		# self._mol = inherit.mol
		# self._db = inherit.db
		
		self._amount = amount
		self._max_lead_time = max_lead_time
		self._supplier = supplier

	### PROPERTIES

	@property
	def db(self):
		return self._db
	
	@property
	def amount(self) -> float:
		"""Returns the amount"""
		return self._amount

	@property
	def id(self):
		return self._compound_id

	@property
	def compound_id(self):
		"""Returns the ID of the associated compound"""
		return self._compound_id

	@property
	def quote_id(self):
		"""Returns the ID of the associated quote"""
		return self._quote_id

	@property
	def max_lead_time(self) -> float:
		"""Returns the max_lead_time from the original quote query"""
		return self._max_lead_time

	@property
	def supplier(self) -> float:
		"""Returns the supplier from the original quote query"""
		return self._supplier

	@amount.setter
	def amount(self, a):
		"""Set the amount and update quotes"""

		quote_id = self.get_cheapest_quote_id( 
			min_amount=a, 
			max_lead_time=self._max_lead_time, 
			supplier=self._supplier,
			none='quiet',
		)

		self._quote_id = quote_id 
		
		self._amount = a

	@property
	def compound(self):
		"""Returns the associated :class:`Compound`"""
		if not self._compound:
			self._compound = self.db.get_compound(id=self.compound_id)
		return self._compound
		
	@property
	def quote(self):
		"""Returns the associated :class:`Quote`"""
		if not self._quote:
			if (q_id := self.quote_id):
				self._quote = self.db.get_quote(id=self.quote_id)
			
			else:
				self._quote = self.compound.get_quotes(
					pick_cheapest=True, 
					min_amount=self.amount, 
					max_lead_time=self.max_lead_time, 
					supplier=self.supplier,
					none='quiet',
				)

		return self._quote

	@property
	def compound_price_amount_str(self):
		# return f'{self} {self.quote.currency_symbol}{self.quote.price} ({self.amount})'
		return f'{self} ({self.amount})'

	@property
	def smiles(self):
		return self.compound.smiles

	@property
	def price(self):
		if self.quote:
			return self.quote.price
		else:
			return None

	@property
	def lead_time(self):
		if self.quote:
			return self.quote.lead_time
		else:
			return None

	### METHODS

	def get_cheapest_quote_id(self, min_amount=None, supplier=None, max_lead_time=None, none='quiet') -> list[dict]:
		"""Query quotes associated to this ingredient"""

		supplier_str = f' AND quote_supplier IS "{supplier}"' if supplier else ""
		lead_time_str = f' AND quote_lead_time <= {max_lead_time}' if max_lead_time else ""
		key_str = f'quote_compound IS {self.compound_id} AND quote_amount >= {min_amount}{supplier_str}{lead_time_str} ORDER BY quote_price'

		result = self.db.select_where(query='quote_id', table='quote', key=key_str, multiple=False, none=none)

		if result:
			quote_id, = result
			return quote_id
		
		else:
			return None

	### DUNDERS

	def __str__(self):
		return f'{self.amount:.2f}mg of C{self._compound_id}'
	
	def __repr__(self):
		return f'{mcol.bold}{mcol.underline}{str(self)}{mcol.unbold}{mcol.ununderline}'

	def __eq__(self, other):

		if self.compound_id != other.compound_id:
			return False

		return self.amount == other.amount
