
import mout
import mcol

from .pose import Pose
from .compound import Compound
from .cset import CompoundSet
from .csetlist import CompoundSetList
from .tset import TagSet
from .pset import PoseSet
from .reaction import Reaction
from .bb import BuildingBlock

from pathlib import Path
import pandas as pd
from pprint import pprint
import json

from rdkit import Chem

from pprint import pprint

from .tools import df_row_to_dict, smiles_has_isotope, clean_smiles

SUPPLIERS = ['enamine', 'mcule']

class HIPPO:

	"""Top-level HIPPO object

	Units
	-----

	* lead time: working days
	* reagant/product quantities: mg
	* pricing: $USD

	SMILES
	------

	HIPPO.Compound's will sanitise their smiles:

		* Stereochemisty removed (@)
		* Only the largest fragment ('.' delimited) is kept

	HIPPO.Pose's will keep stereochemistry in their SMILES fields

	Types of Compound
	-----------------

	Hits:   Experimental structures of protein-ligand complexes.
			The expected format spec is as downloaded from Fragalysis.

			Required inputs:

				- Bound protein-ligand complex PDBs (the ligand is named LIG)

					* The naming convention is:

						TARGET_FRAGMENT_SC_bound.pdb

						TARGET: identifier of the protein target
						FRAGMENT: Five character ID of the fragment e.g. x????
						S: site index (integer)
						C: chain letter (capital)

					* The 'pose_name' is extracted from the name of the pdb 
					  minus the suffix: '_bound.pdb'
				
				- Metadata CSV

					* The relevant row is identified by matching pose_name in 
					  HIPPO to the 'crystal_name' column.
					* The name of the compound is taken from the last five 
					  characters of the 'RealCrystalName' column
					* 'new_smiles' is used if available, else 'smiles'
					* 'RealCrystalName' is stored as Compound.crystal_name
					* 'site_name' is stored as a tag for each Pose object
					* When available 'pdb_entry' becomes Pose.pdb_entry, and
					  the Pose is tagged with 'InPDB'

	Bases/  These are virtual hits, docked/minimised protein-ligand complexes.
	Elabs:  Currently, the expected format is the output of Syndirella.
			
			Required inputs:

				- A minimised .mol file with coordinates for the ligand

				- A minimised bound structure (P+L complex) as a PDB

				- Metadata CSV(s) which are concatenated.
				  They should contain the following information:

					* 'name' which is the unique identifier for a compound
					* 'conformer' which distinguishes different poses / stereoisomers
					* 'smiles' which includes stereochemistry of the conformer
					* 'index' which is shared between conformers of the same compound
					* 'hit_names' space delimited list of fragment inspirations
					* 'ref_pdb' is required only if no bound PDB is provided
					* 'smi_reactant1' for the first reactant (if there is a synthetic route provided)
					* 'smi_reactant2' for the second reactant (if there is a synthetic route provided)

				  Extra metadata if:

					* 'num_atom_difference' the number of atoms added w.r.t. to the base compound
					  == 0 for the base compound itself

	Reactants/          These are not created by the user. The main difference is that they will not 
	Building Blocks:    have poses and their name is the same as the smiles.

	"""

	def __init__(self, project_name, target_name, reactant_catalog='enamine'):

		self._name = project_name
		self._target_name = target_name

		self._max_bb_price = None
		self._min_bb_quantity = None
		self._max_lead_time = None
		
		self._hit_compounds = None

		self._compound_sets = CompoundSetList()
		self._compounds = CompoundSet('compounds')
		# self._poses = PoseSet('poses')
		self._building_blocks = CompoundSet('building_blocks')

		self._protein_feature_strs = None

		assert reactant_catalog in SUPPLIERS
		self._reactant_catalog = reactant_catalog

		self._bbs_id_counter = 0

	### FACTORIES

	@classmethod
	def from_pickle(self, path):
		import pickle
		mout.var('path',str(path),valCol=mcol.file)

		with open(path,'rb') as f:
			self = pickle.load(f)

		return self

	### MAIN METHODS

	def set_cutoffs(self, max_lead_time, max_bb_price, min_bb_quantity):
		self._max_lead_time = max_lead_time
		self._max_bb_price = max_bb_price # do we need this?
		self._min_bb_quantity = min_bb_quantity

	def add_hits(self, metadata_path, pdb_path, pdb_pattern='**/*-x????_??_bound.pdb', tags=None, overwrite=False):
				
		### checks
		
		if not overwrite and 'hits' in [s.name for s in self.compound_sets]:
			mout.error(f'CompoundSet "hits" already exists')
			return

		if not isinstance(metadata_path, Path):
			metadata_path = Path(metadata_path)

		if not isinstance(pdb_path, Path):
			pdb_path = Path(pdb_path)

		tags = tags or ['hits']

		### metadata

		mout.var('metadata_path',str(metadata_path),valCol=mcol.file)

		try:
			metadata_df = pd.read_csv(metadata_path)
		except FileNotFoundError as e:
			mout.error(f'FileNotFoundError: {e}')
			return

		# mout.header('metadata columns:')
		# pprint(list(metadata_df.columns))

		mout.success(f'Parsed metadata CSV')

		### pdbs

		mout.var('pdb_path',str(pdb_path),valCol=mcol.file)
		mout.var('pdb_pattern',pdb_pattern,valCol=mcol.arg)

		pdbs = list(pdb_path.glob(pdb_pattern))

		pdbs = sorted(pdbs)

		if len(pdbs) < 1:
			mout.error(f'Did not find any PDBs',fatal=True,code='HIPPO.add_hits.0')
		
		comp_set = CompoundSet.from_bound_pdbs(name='hits', 
			pdbs=pdbs, 
			metadata_df=metadata_df,
			pdb_pattern=pdb_pattern,
			tags=tags,
			animal=self,
		)

		if 'hits' in self.compound_sets:
			self._compound_sets['hits'] = comp_set
		else:
			self._compound_sets.append(comp_set)
		
		self._hit_compounds = self._compound_sets['hits']

		mout.success(f'Loaded {comp_set.num_compounds} compounds as "{comp_set.name}" ({self.hits.num_poses} poses)')

	def add_elabs(
		self, 
		root_path, 
		tags=None,
		routes_csv_pattern='routes_data/*.csv', 
		elabs_csv_pattern='elabs/*/*/*/*.csv', 
		placements_pattern='elabs/*/*/*/success/*/*.minimised.mol',
		placements_output_dir_name='success',
		minimised_mol_suffix='.minimised.mol',
		elabs_skip_prefix=['.','~$'],
		elabs_skip_substr=['_batch_'],
		elabs_skip_exact=['output.csv'],
		test=False,
		reference_hit=None,
		overwrite=False,
		restart_j=0,
		pickle_dump=None,
		debug=False,
	):

		if not overwrite and 'elabs' in self.compound_sets:
			mout.error(f'CompoundSet "elabs" already exists')
			return
		if not overwrite and 'bases' in self.compound_sets:
			mout.error(f'CompoundSet "bases" already exists')
			return

		if not isinstance(elabs_skip_prefix, list):
			elabs_skip_prefix = [elabs_skip_prefix]
		if not isinstance(elabs_skip_substr, list):
			elabs_skip_substr = [elabs_skip_substr]
		if not isinstance(elabs_skip_exact, list):
			elabs_skip_exact = [elabs_skip_exact]

		tags = tags or ['Syndirella']
		tags = TagSet(tags)
		
		mout.var('root_path',str(root_path))
		mout.var('tags',tags)

		# parse the routes_data CSV
		mout.header('Syndirella Synthetic Routes CSV')
		routes_csv_paths = root_path.glob(routes_csv_pattern)
		mout.var('routes_csv_pattern',routes_csv_pattern)
		mout.var('paths',list(routes_csv_paths))

		# first compound tagged as 'base'?

		# parse the elabs CSV for each batch
		mout.header('Syndirella Elaborations CSVs')
		elabs_csv_paths = list(root_path.glob(elabs_csv_pattern))
		mout.var('elabs_csv_pattern',elabs_csv_pattern)
		mout.var('minimised_mol_suffix',minimised_mol_suffix)
		# mout.var('placements_pattern',placements_pattern)

		# filter
		elabs_csv_paths = [p for p in elabs_csv_paths if not any([p.name.startswith(s) for s in elabs_skip_prefix])]
		elabs_csv_paths = [p for p in elabs_csv_paths if not any([s in p.name for s in elabs_skip_substr])]
		elabs_csv_paths = [p for p in elabs_csv_paths if not any([s == p.name for s in elabs_skip_exact])]

		mout.var('elabs_skip_prefix', elabs_skip_prefix)
		mout.var('elabs_skip_substr', elabs_skip_substr)
		mout.var('elabs_skip_exact', elabs_skip_exact)
		mout.var('#elab CSVs', len(elabs_csv_paths))

		count_bases_without_minimized = 0
		count_bases_missing = 0
		count_isotopic_smiles_skipped = 0
		count_base_with_no_success_or_output = 0
		count_elab_without_base = 0
		count_skipped_multi_name_match = 0

		n_paths = len(elabs_csv_paths)

		# mout.header(f'{n_paths=}')
		# return

		### Create all the Compounds

		# if 'bases' not in self.compound_sets:
			# cset_bases = 

		# # continue from restart (not verified)
		# if restart_j:
		#     if 'bases' in self.compound_sets:
		#         cset_bases = self.compound_sets['bases']
		#     else:
		#         cset_bases = CompoundSet('bases')
		#         self.compound_sets.append(cset_bases)
		#     if 'elabs' in self.compound_sets:
		#         cset = self.compound_sets['elabs']
		#     else:
		#         cset = CompoundSet('elabs')
		#         self.compound_sets.append(cset)

		# # starting afresh
		# else:
		#     cset = CompoundSet('elabs')
		#     cset_bases = CompoundSet('bases')

		#     if overwrite and 'bases' in self.compound_sets:
		#         self.compound_sets['bases'] = cset_bases
		#         mout.debug(f'Overwriting "bases": {cset_bases}')
		#     else:
		#         self.compound_sets.append(cset_bases)
		#         mout.debug(f'Appending "bases": {cset_bases}')

		#     if overwrite and 'elabs' in self.compound_sets:
		#         self.compound_sets['elabs'] = cset
		#         mout.debug(f'Overwriting "elabs": {cset}')
		#     else:
		#         self.compound_sets.append(cset)
		#         mout.debug(f'Appending "elabs": {cset}')

		if debug:
			mout.debug(cset_bases)
			mout.debug(cset)

		for j, path in enumerate(elabs_csv_paths):

			if j < restart_j:
				continue

			# if debug:
			mout.debug(f'{j=}/{n_paths}')

			mout.var('csv',str(path), valCol=mcol.file)

			mout.progress(j, n_paths, append=f'{j=} #elabs={self.num_compounds:>07}')

			df = pd.read_csv(path)

			try:
				names = df['name'].values
			except KeyError:
				mout.error(f'No "name" column in CSV: {path}')
				continue

			if not (path.parent / placements_output_dir_name).is_dir():
				mout.warning(f'Skipping directory with no subdir named "{placements_output_dir_name}" ({path.name})')
				continue

			for i, name in enumerate(names):

				if debug:
					mout.debug(f'{i=}/{len(names)}: {name}')

				if i%100 == 0:
					mout.progress(j, n_paths, append=f'#elabs={self.num_compounds:>07}')

				pose_name = name
				bleach_name = name.replace('_','-')

				# remove conformer info
				name = pose_name[:-2]

				base = 'base' in name

				# comp_output_dir = path.parent / placements_output_dir_name / name
				comp_output_dir = path.parent / placements_output_dir_name / bleach_name

				too_contorted = False

				if not comp_output_dir.is_dir():

					if base:

						if (path.parent / "output" / bleach_name).is_dir():
							comp_output_dir = path.parent / "output" / bleach_name
							if (comp_output_dir / f'{bleach_name}{minimised_mol_suffix}').exists:
								# mout.warning(f'Using too contorted base {pose_name}') 
								count_bases_without_minimized += 1
								too_contorted = True
							else:
								# mout.error(f'Missing base compound {pose_name}')
								count_bases_missing += 1
								continue
						else:
							# mout.error(f'Base has neither "success" nor "output" directory {pose_name}') 
							count_base_with_no_success_or_output += 1
							continue

					else:
						if debug:
							mout.warning(f'No "success" directory {pose_name}') 
							count_elab_without_base += 1
						continue
									
				meta_row = df[df['name'] == pose_name]
				### THIS MAY BE PROBLEMATIC IN THE FUTURE! (missing compounds)
				if len(meta_row) != 0:
					count_skipped_multi_name_match += 1
					continue
				meta_dict = df_row_to_dict(meta_row)

				if smiles_has_isotope(meta_dict["smiles"]):
					count_isotopic_smiles_skipped += 1
					continue

				mol_path = comp_output_dir / f'{bleach_name}{minimised_mol_suffix}'

				if not mol_path.exists:
					continue

				try:
					mol = Chem.MolFromMolFile(str(mol_path))
				except OSError:
					# mout.error(f'Problem w/ mol file for {name}')
					if base:
						count_bases_missing += 1
					continue

				if base:
					comp_tags = tags + ['base']
				else:
					comp_tags = tags + ['elab']

				if too_contorted:
					comp_tags.add('too_contorted')

				if not reference_hit:
					reference_hit = meta_dict['ref_pdb']

				reference_pdb = self.get_hit_pose(reference_hit).pdb_path

				if name not in self.compounds:

					compound = Compound.from_mol(name, mol, tags=comp_tags, animal=self)
					self._add_compound_metadata(compound, meta_dict)

					# if it exists under a different name just add the synthetic route
					if compound in self.compounds:
						if debug:
							mout.debug(f'Adding reaction to elab: {name}')
						self.compounds[compound].add_reaction(compound.reaction)
					else:
						if debug:
							mout.debug(f'Adding elab: {name}')
						self.compounds.add(compound, duplicate='error')
				
				else:
					compound = self.compounds[name]

				if debug:
					mout.debug(f'Adding pose to: {name}')
				
				conformer = meta_dict['conformer']
				site_index = ord(conformer.lower()) - 97
				pose = Pose.from_mol_and_reference(compound, mol, reference_pdb, site_index=site_index, tags=comp_tags)
				pose._stereo_smiles = meta_dict['smiles']
				pose._name = conformer
				compound.add_pose(pose)

			if pickle_dump:
				self.write_pickle(pickle_dump)

			# self.compound_sets['bases'] = cset_bases
			# self.compound_sets['elabs'] = cset

			# mout.debug(cset_bases)
			# mout.debug(cset)

			if test and test == j+1:
				break

		if 'bases' in self.compound_sets:
			self.compound_sets['bases'] = self.get_compounds('base')
		else:
			bases = CompoundSet('bases', self.get_compounds('base'))
			self.compound_sets.append(bases)
			
		if 'elabs' in self.compound_sets:
			self.compound_sets['elabs'] = self.get_compounds('elab')
		else:
			elabs = CompoundSet('elabs', self.get_compounds('elab'))
			self.compound_sets.append(elabs)
		
		self._update_bb_amounts()
		
		mout.finish()

		if pickle_dump:
				self.write_pickle(pickle_dump)

		mout.var('count_bases_without_minimized', count_bases_without_minimized)
		mout.var('count_bases_missing', count_bases_missing)
		mout.var('count_isotopic_smiles_skipped', count_isotopic_smiles_skipped)
		mout.var('count_base_with_no_success_or_output', count_base_with_no_success_or_output)
		mout.var('count_elab_without_base', count_elab_without_base)
		mout.var('count_skipped_multi_name_match', count_skipped_multi_name_match)

		# mout.var('#bases_without_minimized', count_bases_without_minimized)
		# mout.var('#missing_bases', count_bases_missing)
		# mout.var('#isotopic_smiles_skipped', count_isotopic_smiles_skipped)
		# mout.var('#isotopic_smiles_skipped', count_isotopic_smiles_skipped)

		return self.bases, self.elabs

	def quote_reactants(self, quoter):

		from .quoting import NotInCatalogues
		
		count_not_found = 0

		n_bbs = len(self.building_blocks)

		mout.var('#BBs', n_bbs)

		in_cache = False
		last_was_cached = False
		
		# try:
		for i,bb in enumerate(self.building_blocks):

			if bb.price_picker is not None:
				continue
			
			last_was_cached = in_cache
			in_cache = quoter.in_cache(bb)

			if last_was_cached and in_cache:
				if i%50 == 0:
					mout.progress(i, n_bbs, prepend='quoting BBs', append=bb.name)
			else:
				mout.progress(i, n_bbs, prepend='quoting BBs', append=bb.name)
				
			try: 
				bb.tags.add('quote_attempted', duplicate_error=False)
				quoter(bb)
			except NotInCatalogues:
				count_not_found += 1
				continue
			except AssertionError as e:
				raise e
			except KeyboardInterrupt:
				break
			# except Exception as e:
			#     mout.error(f'Quoting error: {e}')
			#     continue

		mout.finish()

		quoter.write_json(f'{self.name}_quoter.json')

		mout.var('#NotInCatalogue',count_not_found)
		mout.var('#Attempted',i)
		mout.success('Finished.')
		
	def add_bases(self, mol_paths, pdb_paths, metadata_paths, tags=None, overwrite=False):

		### checks
		
		if not overwrite and 'bases' in self.compound_sets:
			mout.error(f'CompoundSet "bases" already exists')
			return

		tags = tags or ['bases']
		tags = TagSet(tags)

		comp_set = CompoundSet('bases')

		if isinstance(metadata_paths, list):

			# concatenate all the metadata
			meta_df = None
			for path in metadata_paths:
				df = pd.read_csv(path)

				if meta_df is None:
					meta_df = df
				else:
					meta_df = pd.concat(meta_df, df)

		else:

			meta_df = pd.read_csv(metadata_paths)

		for mol, pdb in zip(mol_paths, pdb_paths):

			if not isinstance(mol, Path):
				mol = Path(mol)

			comp_name = mol.name.replace('.minimised.mol','')

			meta_row = meta_df[meta_df['name'] == comp_name]
			assert len(meta_row) == 1
			meta_dict = df_row_to_dict(meta_row)

			compound = Compound.from_mol(comp_name, mol, tags=tags)

			pose = Pose.from_bound_pdb(compound, pdb)

			compound.add_pose(pose)

			self._add_compound_metadata(compound, meta_dict)

			comp_set.add(compound)

		if overwrite and 'bases' in self.compound_sets:
			self.compound_sets['bases'] = comp_set
		else:
			self.compound_sets.append(comp_set)
		self._base_compounds = self._compound_sets[-1]

		mout.success(f'Loaded {comp_set.num_compounds} compounds as "{comp_set.name}" ({self.bases.num_poses} poses)')

	def random_sample(self, budget, debug=False, max_iter=5000, retry_over_budget=50, increment_count=True):

		# from original HIPPO:
		# by_set = True

		if debug:
			mout.debug(f'animal.random_sample()')

		from .tools import alphabet_index

		# self._bbs_id_counter = 0

		id_str = alphabet_index(self.bbs_id_counter, 3)

		max_lead_time = self.max_lead_time
		mout.var('budget', budget, unit='$USD')
		mout.var('max_lead_time', max_lead_time, unit='days')

		# if debug:
		# 	mout.debug('Getting available reactants...')
		# available_reactants = CompoundSet('available_reactants', [c for c in self.building_blocks if c.price_picker and c.price_picker.min_price], duplicate='no_check')
		# if debug:
		# 	mout.out(available_reactants)

		if debug:
			mout.debug('Getting available products...')
		available_products = CompoundSet('available_products', [c for c in self.compounds if c.reactions], duplicate='no_check')
		if debug:
			mout.debug('Shuffling available products...')
		available_products.shuffle()
		if debug:
			mout.out(available_products)

		bb_set = CompoundSet(f'BBS_{id_str}')
		products = CompoundSet(f'Products of {bb_set.name}')

		if debug:
			mout.header(bb_set)
			mout.header(products)

		iteration = -1
		old_bb_set = None
		over_budget_count = 0

		def cost_of_adding(bb):
			if bb.name in bb_set:
				if debug:
					mout.debug(f'{bb.name} in bb_set')
				now_amount = bb_set[bb.name].amount
				now_cost = bb_set[bb.name].get_price(now_amount)
				delta = bb_set[bb.name].get_price(now_amount + 1) - now_cost
				return delta
			else:
				return bb.get_price(self.min_bb_quantity)

		def add_bb(bb):
			if bb.name in bb_set:
				bb_in_set = bb_set[bb.name]
				bb_in_set.amount += 1
			else:
				bb_set.add(bb.copy(), duplicate='no_check')
				bb_in_set = bb_set[bb.name]
				bb_in_set.amount = self.min_bb_quantity
			if debug:
				mout.var(bb_in_set.name,bb_in_set.amount,unit='mg')

		total_price = 0
		remaining = budget
		while iteration < max_iter:
			
			iteration += 1

			if debug:
				mout.var('iteration', iteration)
				mout.var('total_price', total_price)

			if old_bb_set and over_budget_count > retry_over_budget:
				if total_price > budget - 100:
					mout.success('Good enough')

			if len(available_products) < 1:
				mout.success('Tried all products')
				break

			random_product = available_products.pop()

			if debug:
				mout.debug(f'Trying: {random_product.name}')

			# get cheapest reaction
			reaction_cost_pairs = []
			assert random_product.reactions
			for reaction in random_product.reactions:
				cost = 0
				for bb in reaction.reactants:
					cost += cost_of_adding(bb)
					if not bb.price_picker:
						break
					if not bb.price_picker.min_price:
						break
				else:
					reaction_cost_pairs.append((reaction, cost))

			if len(reaction_cost_pairs) == 0:
				if debug:
					mout.error('no reactions')
				continue

			reaction, adding_cost = sorted(reaction_cost_pairs, key=lambda x: x[1])[0]


			if debug:
				mout.var('cheapest reaction cost', adding_cost, unit='$USD')

			if adding_cost > remaining:
				if debug:
					mout.error('too expensive')
				continue
			elif debug:
				mout.success('cost: OK')

			# add the BBs
			for bb in reaction.reactants:
				add_bb(bb)
			
			# add the product
			product = reaction.product
			product.reaction = reaction
			assert product.num_reactions == 1
			products.add(product)

			# print total BB amount
			if debug:
				mout.var('#BBs', len(bb_set))
				mout.var('#products', len(products))

			# calculate new total cost
			total_price = bb_set.get_price()
			remaining = budget - total_price
			if debug:
				mout.var('total_price', total_price, unit='$USD')
				mout.var('remaining', remaining, unit='$USD')

		else:

			mout.warning(f'Reached {max_iter=}')

			# for bb in random_product.building_blocks:

		mout.var('#products not tried', len(available_products))
		mout.var('total_price', total_price, unit='$USD')
		mout.var('remaining', remaining, unit='$USD')
		mout.var('#BBs', len(bb_set))
		mout.var('#products', len(products))

		if increment_count:
			self._bbs_id_counter += 1

		bb_set.products = products

		return bb_set, products

	def summary(self):
		mout.header(f'{self}')  

		mout.var('target_name', self.target_name)      
		mout.var('max_lead_time', self.max_lead_time, unit='workdays')      
		mout.var('max_bb_price', self.max_bb_price, unit='$')      
		mout.var('min_bb_quantity', self.min_bb_quantity, unit='mg')

		mout.var('#compound_sets',len(self.compound_sets))
		
		all_compounds = self.all_compounds
		mout.var('#compounds',len(all_compounds))
		
		all_poses = self.get_all_poses(all_compounds)
		mout.var('#poses',len(all_poses))
		
		tag_stats = self.get_tag_statistics(all_compounds)

		# mout.var('#building_blocks',len(self.building_blocks))


		if tag_stats:
			mout.var('#tags',len(tag_stats))

		if self.compound_sets:
			mout.out('')
			mout.underline('compound sets:')
			for comp_set in self.compound_sets:
				mout.out(comp_set)
		
			# all_tags = self.all_tags
			if tag_stats:
				mout.out('')
				mout.underline('tags:')
				for tag in tag_stats:

					# num_compounds = len(self.get_compounds(tag))
					# num_poses = len(self.get_poses(tag))

					mout.out(f'{tag} #compounds={tag_stats[tag]["compounds"]}, #poses={tag_stats[tag]["poses"]}')

	def write_pickle(self,file):
		import molparse as mp
		mp.write(file,self)

	### QUERIES

	def get_tag_statistics(self, all_compounds=None):#, all_poses=None):
		all_compounds = all_compounds or self.all_compounds
		# all_poses = all_poses or self.get_all_poses(all_compounds)

		data = {}
		
		for c in all_compounds:

			for t in c.tags:
				if t in data:
					if 'compounds' in data[t]:
						data[t]['compounds'] += 1
					else:
						data[t]['compounds'] = 1
				else:
					data[t] = {}
			
			for p in c.poses:
				for t in p.tags:
					if t in data:
						if 'poses' in data[t]:
							data[t]['poses'] += 1
						else:
							data[t]['poses'] = 1
					else:
						data[t] = {}

		for k,v in data.items():
			if 'compounds' not in v:
				v['compounds'] = 0
			if 'poses' not in v:
				v['poses'] = 0

		return data

	def get_all_poses(self, all_compounds=None):
		all_compounds = all_compounds or self.all_compounds
		poses = sum([[p for p in c.poses] for c in all_compounds],[])
		return PoseSet(poses)

	def get_compounds(self, tag, search=False, inverse=False, all_compounds=None):
		all_compounds = all_compounds or self.all_compounds
		if search:
			return CompoundSet(f'{tag} in tags (search)',[c for c in all_compounds if any([tag in t for t in c.tags])])
		elif inverse:
			assert not search
			return CompoundSet(f'{tag} not in tags',[c for c in all_compounds if tag not in c.tags])
		else:
			return CompoundSet(f'{tag} in tags',[c for c in all_compounds if tag in c.tags])

	def get_poses(self, tag, search=False, inverse=False, all_poses=None, all_compounds=None):
		all_compounds = all_compounds or self.all_compounds
		all_poses = all_poses or self.get_all_poses(all_compounds)
		if search:
			return PoseSet([p for p in all_poses if any([tag in t for t in p.tags])])
		elif inverse:
			assert not search
			return PoseSet([p for p in all_poses if tag not in p.tags])
		else:
			return PoseSet([p for p in all_poses if tag in p.tags])

	def get_hit_pose(self, key):

		""" assume the key is a string of the form: 

			x[0-9][0-9][0-9][0-9]_[0-9][A-Z]

		"""

		comp_name = key.split('_')[0]
		pose_name = key.split('_')[1]

		hit = self.hits[comp_name]

		pose = hit.get_pose(pose_name)

		return pose

	### PLOTTING

	def plot_tag_statistics(self, *args, **kwargs):
		from .plotting import plot_tag_statistics
		return plot_tag_statistics(self, *args, **kwargs)

	def plot_interaction_histogram(self, poses, subtitle=None, **kwargs):
		from .fingerprint import FEATURE_METADATA
		from .plotting import plot_interaction_histogram
		return plot_interaction_histogram(self, poses, FEATURE_METADATA, subtitle, **kwargs)

	def plot_interaction_punchcard(self, poses, subtitle, opacity=1.0, **kwargs):
		from .fingerprint import FEATURE_METADATA
		from .plotting import plot_interaction_punchcard
		return plot_interaction_punchcard(self, poses, FEATURE_METADATA, subtitle, opacity, **kwargs)

	def plot_building_blocks(self, subtitle=None, **kwargs):
		from .plotting import plot_building_blocks        
		return plot_building_blocks(self, subtitle, **kwargs)

	def plot_reactant_amounts(self, subtitle=None, **kwargs):
		from .plotting import plot_reactant_amounts        
		return plot_reactant_amounts(self, subtitle, **kwargs)

	def plot_reactant_price(self, subtitle=None, **kwargs):
		from .plotting import plot_reactant_price        
		return plot_reactant_price(self, subtitle, **kwargs)

	def plot_reactant_sankey(self, subtitle=None, **kwargs):
		from .plotting import plot_reactant_sankey        
		return plot_reactant_sankey(self, subtitle, **kwargs)

	# def plot_reactants_2d(self, subtitle=None, **kwargs):
	#     from .plotting import plot_reactants_2d        
	#     return plot_reactants_2d(self, subtitle, **kwargs)

	def plot_synthetic_routes(self, subtitle=None, cset='elabs', color='num_reactants', **kwargs):
		from .plotting import plot_synthetic_routes        
		return plot_synthetic_routes(self, subtitle, cset, color, **kwargs)

	def plot_numbers(self, subtitle=None, **kwargs):
		from .plotting import plot_numbers        
		return plot_numbers(self, subtitle, **kwargs)

	### INTERNAL METHODS

	def _generate_fingerprints(self, poses):

		n = len(poses)

		fail_count = 0
		for i,pose in enumerate(poses):

			mout.progress(i, n, prepend='fingerprints', append=pose.longname, fill='#')

			try:
				pose.calculate_fingerprint()

			except FailedToAssignBondOrders:
				mout.error(f'Failed to assign bond orders for {pose}')
				pose._fingerprint = None
				fail_count += 1
			# except Exception as e:
			#     mout.error(e)
			#     pose._fingerprint = None
			#     fail_count += 1

		fingerprinted_poses = [c for c in poses if c.fingerprint is not None]

		protein_feature_strs = set(pose.fingerprint.keys()).intersection(*[set(p.fingerprint.keys()) for p in fingerprinted_poses])

		for pose in fingerprinted_poses:
			pose.fingerprint.trim_keys(protein_feature_strs)

		mout.finish()
		if fail_count:
			mout.error(f"Failed to fingerprint {fail_count} poses")
		mout.success(f'Fingerprinted {n - fail_count}/{n} poses')

		mout.var('#fingerprint dimensions',len(protein_feature_strs))

		self._protein_feature_strs = protein_feature_strs

	def _fingerprint_df(self, poses):

		fingerprints = [p.fingerprint for p in poses if p._fingerprint is not None]

		if len(fingerprints) < 1:
			mout.error(f'no fingerprints!')
			return
		
		return pd.DataFrame(fingerprints)
		
	# @mout.debug_log
	def _add_compound_metadata(self, compound, meta_dict, warnings=False):

		### inspirations

		if 'hit_names' in meta_dict:
			for hit_name in meta_dict['hit_names'].split():
				hit_name = hit_name.replace(f'{self.target_name}-','')
				compound.inspirations.append(self.get_hit_pose(hit_name))

		else:
			mout.error(f'No inspiration data for {compound}')

		### base

		if 'Syndirella' in compound.tags and 'base' not in compound.tags:
			
			# strip the conformer
			base_name = meta_dict['base_name']

			# hyphen bleaching
			base_name = base_name.replace('_','-')

			bases = self.get_compounds('base')

			if base_name in bases:
				compound.base = bases[base_name]

		### synthetic routes

		num_reactants = len([k for k in meta_dict if 'reactant' in k and 'smi' in k]) ### OLD FORMAT

		# pprint(meta_dict)

		# mout.debug(num_reactants)

		if num_reactants == 0:
			return

		elif num_reactants == 2:

			# mout.debug('Adding single-step reaction')

			# grab metadata
			reaction_type = meta_dict['reaction']

			if warnings:
				mout.warning('Assuming equal amounts of reactants')
			amounts = 1

			if 'smi_reactant1' in meta_dict:
				reactant1 = meta_dict['smi_reactant1']
				reactant2 = meta_dict['smi_reactant2']
			else:
				reactant1 = meta_dict['reactant1_smi']
				reactant2 = meta_dict['reactant2_smi']

			# create the reactants
			reactant1 = self._get_or_create_building_block(reactant1)
			reactant2 = self._get_or_create_building_block(reactant2)
			
			if reactant1.name_is_smiles:
				metadata_reactant1 = eval(meta_dict['metadata_reactant1'])
				metadata_reactant1 = {d['catalogName']:d for d in metadata_reactant1}
				for v in [v for k,v in metadata_reactant1.items() if self.reactant_catalog in k]:
					reactant1.name = v['catalogId']
					break

			if reactant2.name_is_smiles:
				metadata_reactant2 = eval(meta_dict['metadata_reactant2'])
				metadata_reactant2 = {d['catalogName']:d for d in metadata_reactant2}
				for v in [v for k,v in metadata_reactant2.items() if self.reactant_catalog in k]:
					reactant2.name = v['catalogId']
					break

			# create the reaction
			reaction = Reaction(reaction_type, compound, [reactant1, reactant2], product_amount=amounts, amounts=amounts)

			compound.add_reaction(reaction)

		elif num_reactants > 2:
			mout.error('Not yet implemented')
			raise NotImplementedError

		else:
			mout.error(f'reactants for {compound} = {num_reactants}')
			raise Exception('Unsupported number of reactants')

	def _get_or_create_building_block(self, smiles, debug=False):

		if debug: mout.header(smiles)
		smiles = clean_smiles(smiles, verbosity=debug)['smiles']

		if smiles in self.building_blocks.smiles:
			if debug: mout.debug('returning existing BB')
			return self.building_blocks.get_by_smiles(smiles)

		else:
			bb = BuildingBlock(smiles)
			self.building_blocks.add(bb)
			if debug: mout.debug('returning new BB')
			return bb

	def _get_or_create_compound(self, name, smiles, debug=False, duplicate='error'):

		if debug: mout.header(name)
		# smiles = clean_smiles(smiles, verbosity=debug)['smiles']

		if name in self.compounds.names:
			if debug: mout.debug('returning existing Compound')
			return self.compounds.get_by_name(name)

		else:
			comp = Compound(name, smiles)
			self.compounds.add(comp, duplicate=duplicate)
			if debug: mout.debug('returning new Compound')
			return comp

	# def _get_or_create_pose(self, longname, debug=False):

	#     if debug: mout.header(name)
	#     # smiles = clean_smiles(smiles, verbosity=debug)['smiles']

	#     if name in self.compounds.names:
	#         if debug: mout.debug('returning existing Compound')
	#         return self.compounds.get_by_name(name)

	#     else:
	#         comp = Compound(name, smiles)
	#         self.compounds.add(comp)
	#         if debug: mout.debug('returning new Compound')
	#         return comp

	def _update_bb_amounts(self, debug=False):

		for bb in self.building_blocks:
			bb.amount = 0

		# if debug:
			# mout.debug('getting Syndirella compounds')
		# syndirella_comps = self.get_compounds('Syndirella')

		n = len(self.compounds)
		mout.var('#comps', n)
		
		for i, comp in enumerate(self.compounds):
			if debug and i%100 == 0:
				mout.progress(i,n)
			for reaction in comp.reactions:
				for reactant in reaction.reactants:
					self.building_blocks[reactant].amount += 1
		if debug:
			mout.finish()

	def _delete_compounds(self, delete_list, debug=False):

		delete_references = 0
		for cset in self.compound_sets:
			if debug:
				mout.debug(cset)
				n = len(cset)
			for i,comp in enumerate(cset):
				if debug and i%100 == 0:
					mout.progress(i,n)
				if comp in delete_list:
					del cset[comp]
					delete_references += 1
			if debug:
				mout.finish()
					
		delete_compounds = 0
		n = self.num_compounds
		if debug:
			mout.debug('Deleting from animal.compounds')
			mout.debug(self.compounds)
		for i,comp in enumerate(self.compounds):
			if debug and i%100 == 0:
				mout.progress(i,n)
			if comp in delete_list:
				del self.compounds[comp]
				delete_compounds += 1

		if debug:
			mout.finish()
			mout.debug(self.compounds)

		mout.var('references deleted', delete_references)
		mout.var('compounds deleted', delete_compounds)

	def _delete_reactants(self, delete_list, debug=False):
					
		delete_compounds = 0
		n = self.num_compounds
		if debug:
			mout.debug('Deleting from animal.building_blocks')
			mout.debug(self.building_blocks)
		for i,comp in enumerate(self.building_blocks):
			if debug and i%100 == 0:
				mout.progress(i,n)
			if comp in delete_list:
				del self.building_blocks[comp]
				delete_compounds += 1

		if debug:
			mout.finish()
			mout.debug(self.building_blocks)

		mout.var('reactants deleted', delete_compounds)

	### PROPERTIES

	@property
	def compounds(self):
		return self._compounds

	@property
	def name(self):
		return self._name

	@property
	def num_compounds(self):
		return len(self.all_compounds)

	@property
	def all_compounds(self):
		# if len(self.compound_sets) == 1:
			# return self.compound_sets[0]
		# return set().union(*self.compound_sets)
		# return self.compound_sets.all_compounds
		return self.compounds

	@property
	def num_poses(self):
		return len(self.all_poses)

	@property
	def all_poses(self):
		return self.get_all_poses()
		# mout.debug('v2')

		# return PoseSet(sum([c.poses for c in self.all_compounds], PoseSet()))

	@property
	def all_tags(self):

		all_tags = TagSet()

		for thing in self.all_compounds:
			for tag in thing.tags:
				if tag not in all_tags:
					all_tags.add(tag)

		for thing in self.all_poses:
			for tag in thing.tags:
				if tag not in all_tags:
					all_tags.add(tag)

		return all_tags

	@property
	def num_tags(self):
		return len(self.all_tags)
	
	@property
	def compound_sets(self):
		return self._compound_sets

	@property
	def max_lead_time(self):
		return self._max_lead_time

	@property
	def max_bb_price(self):
		return self._max_bb_price

	@property
	def min_bb_quantity(self):
		return self._min_bb_quantity

	@property
	def hits(self):
		return self.compound_sets['hits']

	@property
	def bases(self):
		if 'bases' in self.compound_sets:
			return self.compound_sets['bases']
		else:
			return self.get_compounds('base')
	
	@property
	def target_name(self):
		return self._target_name

	@property
	def building_blocks(self):
		return self._building_blocks

	@property
	def num_building_blocks(self):
		return len(self.building_blocks)

	@property
	def elabs(self):
		if 'elabs' in self.compound_sets:
			return self.compound_sets['elabs']
		else:
			return self.get_compounds('elab')

	@property
	def reactant_catalog(self):
		return self._reactant_catalog

	@property
	def bbs_id_counter(self):
		return self._bbs_id_counter
	
	### DUNDERS

	def __repr__(self):
		return f'HIPPO({self.name})'

class FailedToAssignBondOrders(Exception):
	pass
