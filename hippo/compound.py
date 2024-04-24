
import mcol

from rdkit import Chem

from .pose import Pose
# from .pset import PoseSet
from .tags import TagSet
# from .rset import ReactionSet
from .target import Target

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
		self._num_heavy_atoms = None

		self._metadata = metadata
		
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

		serialisable_fields = ['id','inchikey','alias','smiles']

		data = {}
		for key in serialisable_fields:
			data[key] = getattr(self, key)

		if base := self.base:
			data['base'] = base.inchikey

		if metadata := self.metadata:
			for key in metadata:
				data[key] = metadata[key]

		return data

	@property
	def table(self):
		return self._table
	
	
	### METHODS

	def get_tags(self) -> set:
		tags = self.db.select_where(query='tag_name', table='tag', key='compound', value=self.id, multiple=True, none='quiet')
		return TagSet(self, {t[0] for t in tags}, commit=False)

	def get_quotes(self, min_amount=None, supplier=None, max_lead_time=None, none='quiet', pick_cheapest=False, df=False) -> list[dict]:
		"""Get all quotes associated to this compound"""

		quote_ids = self.db.select_where(query='quote_id', table='quote', key='compound', value=self.id, multiple=True, none=none)

		if quote_ids:
			quotes = [self.db.get_quote(id=q[0]) for q in quote_ids]
		else:
			return []

		if supplier:
			quotes = [q for q in quotes if q.supplier == supplier]

		if min_amount:
			quotes = [q for q in quotes if q.amount >= min_amount]

		if max_lead_time:
			quotes = [q for q in quotes if q.lead_time <= max_lead_time]

		if pick_cheapest:
			return sorted(quotes, key=lambda x: x.price)[0]

		if df:
			from pandas import DataFrame
			return DataFrame([q.asdict() for q in quotes]).drop(columns='compound')
		
		return quotes

	def get_reactions(self, none='error', as_reactant=False) -> list:
		"""Get the associated reactions as product, unless as_reactant is True."""

		from .rset import ReactionSet

		if as_reactant:
			reaction_ids = self.db.select_where(query='reactant_reaction', table='reactant', key='compound', value=self.id, multiple=True, none=none)
		else:
			reaction_ids = self.db.select_where(query='reaction_id', table='reaction', key='product', value=self.id, multiple=True, none=none)

		reaction_ids = [q for q, in reaction_ids]

		return ReactionSet(self.db, reaction_ids)

	def get_poses(self, 
		target: str = None
	) -> list[Pose]:

		pose_ids = self.db.select_where(query='pose_id', table='pose', key='compound', value=self.id, multiple=True, none=False)

		if not pose_ids:
			return None

		# poses = [self.db.get_pose(id=q[0]) for q in pose_ids]

		from .pset import PoseSet

		return PoseSet(self.db, [q[0] for q in pose_ids])

	def get_dict(self, mol=True, reactions=False, metadata=True, count_by_target=False):
		
		"""Returns a dictionary representing this Compound"""

		serialisable_fields = ['id','alias', 'inchikey', 'smiles', 'num_poses', 'num_reactant', 'num_reactions']

		# poses
		# reactions
		# poses.targets

		assert not reactions

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
		
		poses = self.poses

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
		else:
			quote = quote.id
		
		return Ingredient(self.db, self.id, amount, quote, supplier, max_lead_time)

	def draw(self, align_substructure: bool = False):
		"""Display this compound (and its base if it has one)"""

		if self.base:
			from molparse.rdkit import draw_mcs
			drawing = draw_mcs({self.base.smiles:f'{self.base} (base)', self.smiles:str(self)}, align_substructure=align_substructure, show_mcs=False, highlight_common=False)
			display(drawing)
		else:
			display(self.mol)

	def summary(self, metadata: bool = True, draw: bool = True,):
		"""Print a summary of this compound"""
		logger.header(repr(self))
		logger.var('inchikey', self.inchikey)
		logger.var('alias', self.alias)
		logger.var('smiles', self.smiles)
		logger.var('base', self.base)
		poses = self.poses
		logger.var('#poses', len(poses))
		if poses:
			logger.var('targets', poses.targets)
		logger.var('#reactions (product)', self.num_reactions)
		logger.var('#reactions (reactant)', self.num_reactant)
		logger.var('tags', self.tags)
		if draw:
			self.draw()
		if metadata:
			logger.var('metadata', str(self.metadata))

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
    
		if isinstance(compound, int):
			self._compound_id = compound
		else:
			self._compound_id = compound.id
			self._compound = None

		if isinstance(quote, int):
			self._quote_id = quote
			self._quote = None
		elif quote is None:
			self._quote_id = None
			self._quote = None
		else:
			self._quote_id = quote.id
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

		quote = self.get_quotes(
			pick_cheapest=True, 
			min_amount=a, 
			max_lead_time=self._max_lead_time, 
			supplier=self._supplier,
			none='quiet',
		)
		
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
		if not self._quote and (q_id := self.quote_id):
			self._quote = self.db.get_quote(id=self.quote_id)
		return self._quote

	### DUNDERS

	def __str__(self):
		return f'{self.amount:.2f}mg of C{self._compound_id}'
	
	def __repr__(self):
		return f'{mcol.bold}{mcol.underline}{str(self)}{mcol.unbold}{mcol.ununderline}'

	def __eq__(self, other):

		if self.compound_id != other.compound_id:
			return False

		return self.amount == other.amount
