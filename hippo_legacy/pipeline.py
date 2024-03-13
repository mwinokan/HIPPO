
import mout

mout.debug('imports...')
import mcol
import molparse as mp

import os

# silence NumbaDeprecationWarning
import warnings
from numba.core.errors import NumbaDeprecationWarning
warnings.filterwarnings("ignore", category=NumbaDeprecationWarning)

import pandas as pd

from .set import CompoundSet
from .compound import Compound
from .block import BuildingBlock, BuildingBlockSet

import pickle

from pprint import pprint

import plotly.graph_objects as go
import plotly.express as px

from .compound import FailedToAssignBondOrders
import numpy as np

import json
import random

import glob

# from scipy.special import seterr
# seterr(all='raise')

from numpy import seterr
seterr(divide='raise')

class HIPPO:

    """🦛 Hit Interaction Profiling for Procurement Optimisation

    Pipeline:
    ---------

    1. Add protein reference PDB: 
        add_protein_reference()

    2. Add crystallographic data (aligned folder): 
        add_hit_directory()
    
    3. Add compound data:
        add_compound_directory() 

    """

    ### DUNDERS

    def __init__(self, project_key, verbosity=2):
        self._project_key = project_key
        self._verbosity = verbosity

        self._protein_reference_path = None
        self._protein_system = None
        self._hit_pdbs = None
        self._metadata_df = None
        self._compound_sets = []
        self._missed_features = None
        self._building_blocks = BuildingBlockSet()
        self._hit_features = None
        self._bb_id_to_smiles = None
        self._bb_smiles_to_id = None
        self._hit_compounds = None
        self._base_compounds = None

        # stage flags
        self._has_protein_reference = False
        self._has_hit_metadata = False
        self._has_hit_directory = False
        self._has_all_compounds = False
        self._has_exported_bb_smiles = False
        self._has_purchase_info = False
        # self._has_been_culled_by_lead_time = False
        self._has_products = False
        self._has_hit_fingerprints = False
        self._has_product_fingerprints = False
        self._has_non_purchaseable_been_culled = False
        self._has_bb_map = None

        # umap
        self._umap_dims = None
        self._reducer = None
        self._scaler = None
        self._transformer = None
        self._embedding = None
        self._umap_df = None

        self._comp_df_meta_keys = [
            'name',
            'smiles',
            'set_name',
            'is_pains',
            'cost_range_str',
            'lead_time',
        ]

    def __repr__(self):
        return f'HIPPO({self.project_key})'

    ### FACTORIES

    @classmethod
    def from_pickle(self, path, verbosity=3):
        if verbosity > 1:
            mout.debug('HIPPO.from_pickle()')
            mout.var('path',str(path),valCol=mcol.file)

        with open(path,'rb') as f:
            self = pickle.load(f)

        self._verbosity = verbosity

        return self

    ### PROPERTIES

    @property
    def project_key(self):
        return self._project_key

    @property
    def verbosity(self):
        return self._verbosity

    @property
    def protein_system(self):
        return self._protein_system

    @property
    def metadata_df(self):
        return self._metadata_df

    @property
    def compound_sets(self):
        return self._compound_sets

    @property
    def all_compounds(self):
        return set().union(*self.compound_sets)

    @property
    def fingerprinted_compounds(self):
        return [c for c in self.all_compounds if c.fingerprint is not None]

    @property
    def num_compounds(self):
        return len(self.all_compounds)

    @property
    def num_products(self):
        return len(self.all_products)

    @property
    def num_building_blocks(self):
        return len(self.building_blocks)

    @property
    def building_blocks(self):
        return self._building_blocks

    @property
    def compound_df(self):
        
        data = []

        for comp in self.fingerprinted_compounds:

            data_row = dict(comp.fingerprint)

            data.append(comp.dict_with_fingerprint)

        return pd.DataFrame(data)

    @property
    def has_protein_reference(self):
        return self._has_protein_reference

    @property
    def has_hit_metadata(self):
        return self._has_hit_metadata

    @property
    def has_hit_directory(self):
        return self._has_hit_directory

    @property
    def has_all_compounds(self):
        return self._has_all_compounds
    
    @has_all_compounds.setter
    def has_all_compounds(self, b):
        self._has_all_compounds = b

    @property
    def has_exported_bb_smiles(self):
        return self._has_exported_bb_smiles
    
    @has_exported_bb_smiles.setter
    def has_exported_bb_smiles(self, b):
        self._has_exported_bb_smiles = b

    # @property
    # def has_been_culled_by_lead_time(self):
        # return self._has_been_culled_by_lead_time

    @property
    def has_products(self):
        return self._has_products

    @property
    def has_hit_fingerprints(self):
        return self._has_hit_fingerprints

    @property
    def has_product_fingerprints(self):
        return self._has_product_fingerprints

    @property
    def has_non_purchaseable_been_culled(self):
        return self._has_non_purchaseable_been_culled

    @property
    def has_purchase_info(self):
        return self._has_purchase_info

    @has_purchase_info.setter
    def has_purchase_info(self,b):
        self._has_purchase_info = b

    @property
    def hit_features(self):
        if self._hit_features is None:
            self._hit_features = self.compound_sets[0].get_present_features()
        return self._hit_features

    @property
    def has_bb_map(self):
        return self._has_bb_map

    @property
    def hits(self):
        return self._hit_compounds

    @property
    def bases(self):
        return self._base_compounds

    @property
    def product_sets(self):
        assert self.hits is not None
        assert self.has_hit_directory
        return self._compound_sets[1:]

    ### METHODS

    def write_pickle(self,file):
        mp.write(file,self)
    
    def add_protein_reference(self, path):

        if self.has_protein_reference:
            mout.error(f'{self} already has a protein reference')
            return

        if self.verbosity > 1:
            mout.debug('HIPPO.add_protein_reference()')
            mout.var('path',str(path),valCol=mcol.file)

        self._protein_reference_path = path
        self._protein_system = mp.parsePDB(path,verbosity=self.verbosity-2,alternative_site_warnings=False)

        self.protein_system.name = f'{self.project_key} [protein reference]'

        if self.verbosity > 1:
            self.protein_system.summary()

        for chain in self.protein_system.chains:
            assert chain.type == 'PRO'

        if self.verbosity:
            mout.success(f'Parsed protein PDB')

        self._has_protein_reference = True

    def add_hit_metadata(self, path):

        assert self.has_protein_reference

        if self.has_hit_metadata:
            mout.error(f'{self} already has hit metadata')
            return

        if self.verbosity > 1:
            mout.debug('HIPPO.add_hit_metadata()')
            mout.var('path',str(path),valCol=mcol.file)

        self._metadata_df = pd.read_csv(path)

        if self.verbosity > 2:
            mout.header('metadata columns:')
            pprint(list(self.metadata_df.columns))

        if self.verbosity:
            mout.success(f'Parsed metadata CSV')

        self._has_hit_metadata = True

    def add_hit_directory(self, name='hits', path=None, pattern='*-x????_??_bound.pdb'):

        assert self.has_hit_metadata

        if self.has_hit_directory:
            mout.error(f'{self} already has hit directory')
            return

        assert path is not None

        if name in [s.name for s in self.compound_sets]:
            mout.warning(f'Skipping existing {name}')

        if self.verbosity > 1:
            mout.debug('HIPPO.add_hit_data()')
            mout.var('path',str(path),valCol=mcol.file)
            mout.var('pattern',pattern,valCol=mcol.arg)

        pdbs = list(path.glob(f'**/{pattern}'))

        pdbs = sorted(pdbs)

        if len(pdbs) < 1:
            mout.error(f'Did not find any PDBs',fatal=True,code='HIPPO.add_hit_data.0')

        self._hit_pdbs = pdbs
        
        comp_set = CompoundSet.from_bound_pdbs(name='hits', 
            pdbs=self._hit_pdbs, 
            metadata_df=self.metadata_df,
            prefix=f'{self.project_key}-'
        )

        self._compound_sets.append(comp_set)
        self._hit_compounds = self._compound_sets[-1]

        if self.verbosity:
            mout.success(f'Loaded {comp_set.num_compounds} compounds "{comp_set.name}"')

        self._has_hit_directory = True

    def add_compounds_from_sdf(self, name, path, idName='name',molColName='mol'):

        assert self.has_hit_directory

        if self.has_all_compounds:
            mout.error(f'{self} is tagged as having all compounds')
            return

        if name in [s.name for s in self.compound_sets]:
            mout.warning(f'Skipping existing {name}')

        if self.verbosity > 1:
            mout.debug('HIPPO.add_compounds_from_sdf()')
            mout.var('name',name,valCol=mcol.arg)
            mout.var('path',str(path),valCol=mcol.file)
            mout.var('idName',idName,valCol=mcol.arg)
            mout.var('molColName',molColName,valCol=mcol.arg)
        
        from rdkit.Chem import PandasTools
        mol_df = PandasTools.LoadSDF(str(path), idName=idName, molColName=molColName)

        if self.verbosity > 2:
            mout.header('metadata columns:')
            pprint(list(mol_df.columns))

        comp_set = CompoundSet.from_df(name=name, df=mol_df[[idName,molColName]], protein=self.protein_system, verbosity=self.verbosity-1)
        
        if self.verbosity > 3:
            print(comp_set)
            pprint(comp_set.compounds)

        self._compound_sets.append(comp_set)

        if self.verbosity:
            mout.success(f'Loaded {comp_set.num_compounds} compounds "{comp_set.name}"')

    def add_reaction_products(self, name, metadata, data_path, mol_pattern, skip_unconstrained_minimisation_fail=False, expecting=None, skip_too_moved=True, skip_too_contorted=True):

        assert self.has_hit_directory

        if self.has_all_compounds:
            mout.error(f'{self} is tagged as having all compounds')
            return

        if name in [s.name for s in self.compound_sets]:
            mout.warning(f'Skipping existing {name}')

        if self.verbosity > 1:
            mout.debug('HIPPO.add_compound_set()')
            mout.var('name',name,valCol=mcol.arg)
            mout.var('metadata',str(metadata),valCol=mcol.file)
            mout.var('data_path',str(data_path),valCol=mcol.file)
            mout.var('mol_pattern',mol_pattern,valCol=mcol.arg)

        meta_df = pd.read_csv(metadata)
            
        if self.verbosity > 2:
            mout.header('metadata columns:')
            pprint(list(meta_df.columns))

        ### get the building blocks

        for index, vals in meta_df[['reactant1_smi','reactant1_mw','reactant1_metadata','reactant2_smi','reactant2_mw','reactant2_metadata']].iterrows():

            smiles1, mol_weight1, metadata1, smiles2, mol_weight2, metadata2 = vals

            # skip empty rows
            if not isinstance(smiles1, str):
                continue

            if not isinstance(smiles2, str):
                continue

            if smiles1 not in self.building_blocks:

                ### REACTANT 1

                try:
                    mol_weight1 = float(mol_weight1)
                except ValueError:
                    mout.error(f"{index=} Can't convert str-->float (mol_weight1)")

                bb = BuildingBlock(smiles1, mol_weight1) #, metadata1)

                if bb not in self.building_blocks:
                    self.building_blocks.add(bb)

            if smiles2 not in self.building_blocks:
                
                ### REACTANT 2
            
                try:
                    mol_weight2 = float(mol_weight2)
                except ValueError:
                    mout.error(f"{index=} Can't convert str-->float (mol_weight2)")

                bb = BuildingBlock(smiles2, mol_weight2) #, metadata2)

                if bb not in self.building_blocks:
                    self.building_blocks.add(bb)

        assert self.building_blocks

        if self.verbosity:
            mout.success(f'Found {self.num_building_blocks} building blocks for "{name}"')

        comp_set = CompoundSet(name=name)

        mols = list(data_path.glob(mol_pattern))

        basename_mols = [os.path.basename(m) for m in mols]

        if expecting and len(mols) != expecting:
            mout.error(f'{name}: {len(mols)=} != {expecting}')

            if 'name' in meta_df:

                for meta_name in meta_df['name'].values:
                    if not isinstance(meta_name,str):
                        continue

                    if f'{meta_name}.minimised.mol' not in basename_mols:
                        mout.warning(f'Missing {meta_name}.minimised.mol')

        # get the base compound
        base = list(data_path.glob('base/base.minimised.mol'))
        assert base

        progress = len(mols) > 100

        for i,mol_file in enumerate(base + mols):

            prod_name = mol_file.name.split('.')[0]

            if progress:
                mout.progress(i,len(mols),prepend=name,append=prod_name,fill='#')

            mol = mp.parse(mol_file, verbosity=0)

            compound = Compound.from_rdkit_mol(prod_name, mol)

            compound._metadata_csv = str(metadata)

            if prod_name == 'base':
                place_type = 'frag'
            elif 'base' in prod_name:
                place_type = 'base'
            elif 'frag':
                place_type = 'frag'
            else:
                raise Exception

            if f'{place_type}_placed_name' in meta_df:
                meta_row = meta_df[meta_df[f'{place_type}_placed_name'] == prod_name]
                assert len(meta_row), prod_name
            else:
                meta_row = meta_df[meta_df['name'] == prod_name]
                assert len(meta_row), prod_name
            
            if f'{place_type}_fail_minimization_without_constraint' in meta_row:
                compound._minimization_failed_without_constraint = meta_row[f'{place_type}_fail_minimization_without_constraint'].values[0]

            if 'product_smi' in meta_row:
                compound._smiles = str(meta_row['product_smi'].values[0])
            elif 'smiles' in meta_row:
                compound._smiles = str(meta_row['smiles'].values[0])
            else:
                mout.error(f"Couldn't get SMILES from metadata {index=}")
                raise Exception("No SMILES")

            if 'is_pains' in meta_row:
                compound._is_pains = bool(meta_row['is_pains'].values[0])
            else:
                mout.error(f"Couldn't get PAINS from metadata {index=}")
                raise Exception("No PAINS")

            if 'num_atom_difference' in meta_row:
                compound._num_atom_difference = int(meta_row['num_atom_difference'].values[0])

            if '∆∆G' in meta_row:
                compound._fragmenstein_ddG = float(meta_row['∆∆G'].values[0])
            else:
                # get fragmenstein json
                try:
                    with open(str(mol_file).replace('.mol','.json')) as f:
                        mol_json_data = json.load(f)
                        # mout.json(mol_json_data)

                    compound._fragmenstein_ddG = float(mol_json_data['Energy']['ligand_ref2015']['total_score'] - mol_json_data['Energy']['unbound_ref2015']['total_score'])
                    compound._fragmenstein_mRMSD = float(mol_json_data['mRMSD'])
                    compound._fragmenstein_too_moved = compound._fragmenstein_mRMSD > 1.0
                    compound._fragmenstein_too_contorted = compound._fragmenstein_ddG > 0.0

                    if skip_too_moved and compound._fragmenstein_too_moved:
                        # mout.warning(f'Skipping {compound} [too moved]')
                        continue

                    if skip_too_contorted and compound._fragmenstein_too_contorted:
                        # mout.warning(f'Skipping {compound} [too contorted]')
                        continue

                except FileNotFoundError:
                    pass

            if 'LE' in meta_row:
                compound._fragmenstein_ligand_efficiency = meta_row['LE'].values[0]

            if 'outcome' in meta_row:
                compound._fragmenstein_outcome = str(meta_row['outcome'].values[0])

            compound._protein_system = self.protein_system
            
            reactant1_smiles = str(meta_row['reactant1_smi'].values[0])
            reactant2_smiles = str(meta_row['reactant2_smi'].values[0])
            
            bb1 = self.building_blocks[reactant1_smiles]
            bb2 = self.building_blocks[reactant2_smiles]
            compound._building_blocks = BuildingBlockSet([bb1,bb2])

            if 'base' == prod_name:
                if self.bases is None:
                    assert len(self.compound_sets) == 1
                    self._base_compounds = CompoundSet('bases')
                    self._compound_sets.append(self.bases)

                compound._name = f'{name}-base'
                self.bases.add(compound)

            else:
                comp_set.add(compound)
        
        mout.finish()

        self._compound_sets.append(comp_set)

        if self.verbosity:
            mout.success(f'Loaded {comp_set.num_compounds} compounds "{comp_set.name}"')

    def generate_fingerprints(self,compounds):

        if self.verbosity > 1:
            mout.debug('HIPPO.generate_fingerprints()')

        n = len(compounds)

        fail_count = 0
        for i,compound in enumerate(compounds):
            if self.verbosity:
                mout.progress(i, n, prepend='fingerprints', append=f'{mcol.bold}{compound.set_name}{mcol.clear} {compound.name}', fill='#')

            try:
                compound.calculate_fingerprint()

            except FailedToAssignBondOrders:
                mout.error(f'Failed to assign bond orders for {compound}')
                compound._fingerprint = None
                fail_count += 1
            except Exception as e:
                mout.error(e)
                compound._fingerprint = None
                fail_count += 1

        fingerprinted_compounds = [c for c in compounds if c.fingerprint if not None]

        protein_feature_strs = set(compound.fingerprint.keys()).intersection(*[set(c.fingerprint.keys()) for c in fingerprinted_compounds])

        for compound in fingerprinted_compounds:
            compound.fingerprint.trim_keys(protein_feature_strs)

        if self.verbosity > 1:
            mout.var('#feature_dims',len(protein_feature_strs))

        if self.verbosity:
            mout.finish()
            if fail_count:
                mout.error(f"Failed to fingerprint {fail_count} compounds")
            mout.success(f'Fingerprinted {n - fail_count}/{n} compounds')

        return protein_feature_strs

    def score_interaction_coverage(self, show=False):

        if self.verbosity > 1:
            mout.debug('HIPPO.score_interaction_coverage()')

        self._missed_features = {}

        features_by_set = {}
        for comp_set in self.compound_sets:
            features_by_set[comp_set.name] = comp_set.get_present_features()

        for i,comp_set1 in enumerate(self.compound_sets):
            for comp_set2 in self.compound_sets[i+1:]:

                name1 = comp_set1.name
                name2 = comp_set2.name
            
                features1 = features_by_set[name1]
                features2 = features_by_set[name2]

                shared = features1.intersection(features2)

                mout.header(f'intersection({name1}, {name2})')
                
                in_1_but_not_2, in_2_but_not_1 = self.compare_feature_sets(name1,name2,features1,features2)

                num_shared = len(shared)
                num_in_1_but_not_2 = len(in_1_but_not_2)
                num_in_2_but_not_1 = len(in_2_but_not_1)

                self._missed_features[f'in {name1} but not {name2}'] = in_1_but_not_2
                self._missed_features[f'in {name2} but not {name1}'] = in_2_but_not_1

                ### SANKEY PLOT

                fig = go.Figure()

                trace = go.Sankey(
                    node = dict(
                        pad = 15,
                        thickness = 20,
                        # line = dict(color = "black", width = 0.5),
                        label = [name1, f"only in {name2}", f"only in {name1}", f"shared"],
                        # color = "blue"
                    ),
                    link = dict(
                        source = [0, 0, 1], # indices correspond to labels, eg A1, A2, A1, B1, ...
                        target = [2, 3, 3],
                        value =  [num_in_1_but_not_2, num_shared, num_in_2_but_not_1]
                    )
                )

                fig.add_trace(trace)

                if show:
                    fig.show()

                # write the figure
                mp.write(f'sankey_{name1}_{name2}.html',fig)

    def compare_feature_sets(self, name1, name2, set1, set2, print_diff=False, verbosity=1):

        assert set1
        assert set2

        in_1_but_not_2 = set1 - set2
        if verbosity:
            mout.var(f'# only in {name1}',len(in_1_but_not_2))
            mout.var(f'% missed {name1} features',f'{len(in_1_but_not_2)/len(set1):.1%}')
            mout.var(f'% covered {name1} features',f'{(len(set1)-len(in_1_but_not_2))/len(set1):.1%}')
            if print_diff:
                pprint(in_1_but_not_2)

        in_2_but_not_1 = set2 - set1
        if verbosity: 
            mout.var(f'# only in {name2}',len(in_2_but_not_1))
            mout.var(f'% missed {name2} features',f'{len(in_2_but_not_1)/len(set2):.1%}')
            mout.var(f'% covered {name2} features',f'{(len(set2)-len(in_2_but_not_1))/len(set2):.1%}')
            if print_diff:
                pprint(in_2_but_not_1)

        return in_1_but_not_2, in_2_but_not_1

    def random_building_block_set(self, budget=None, size=None, max_lead_time=None, verbosity=1, max_iter=5000, start_with=None, top_down=False, debug=False, by_set=False, bbs_id=None, bb_price_cutoff=None, retry_over_budget=50):

        assert size or budget, "One of size and budget must be specified"
        assert not (size and budget), "One of size and budget must be None"

        # if max_lead_time:
        #     available = [bb for bb in self.building_blocks if bb.lead_time <= max_lead_time]
        #     if verbosity:
        #         mout.warning(f'Ignoring {1 - len(available)/len(self.building):.1%} with lead times exceeding {max_lead_time} weeks')
        # else:
        available = self.building_blocks

        if size:

            mout.out(f'Choosing {size} from {len(available)} building blocks')

            # get a random sample of building blocks
            bb_set = BuildingBlockSet(random.sample(self.building_blocks,size))

            bb_set.name = 'random'

        elif budget:

            if not top_down:
            
                if max_lead_time:
                    mout.out(f'Choosing from {len(available)} building blocks up to ${budget} with lead_time <= {max_lead_time} days')
                else:
                    mout.out(f'Choosing from {len(available)} building blocks up to ${budget}')

                if start_with:

                    smiles = random.sample(start_with,1)
                    bb = self.building_blocks[smiles]
                    bb_set = BuildingBlockSet(bb)

                else:
                    # start with a random BB
                    bb_set = BuildingBlockSet(random.sample(self.building_blocks,1))
                
                mout.var('start',bb_set.smiles)

                all_products = self.all_products
                all_compounds = self.all_compounds

                iteration = 0

                while iteration < max_iter:

                    if verbosity > 1:
                        mout.var(f'{mcol.bold}iteration',iteration)
                        mout.var('#bbs',len(bb_set))

                    # add a BB
                    available = self.building_blocks - bb_set
                    addition = random.sample(available,1)[0]
                    if verbosity > 1:
                        mout.out(f'Adding: {addition.smiles}')
                    bb_set.add(addition)

                    # get products
                    products = bb_set.get_products(all_compounds)
                    if verbosity > 1:
                        mout.var('#products',len(products))

                    # calculate total cost & lead time
                    total_price = bb_set.get_price()
                    remaining = budget - total_price
                    if verbosity > 1:
                        mout.var('total_price',total_price)
                        mout.var('budget remaining',remaining)

                    lead_times = bb_set.get_lead_times()
                    this_max_lead_time = max(lead_times)
                    if verbosity > 1:
                        mout.var('avg(lead_times)',np.mean(lead_times))
                        if max_lead_time:
                            mout.var('max(lead_times)',this_max_lead_time,valCol=mcol.error if this_max_lead_time > max_lead_time else None)
                        else:
                            mout.var('max(lead_times)',this_max_lead_time)

                    if remaining < 0:
                        if verbosity > 1:
                            mout.error('Over budget')
                        bb_set = old_bb_set.copy()
                        break

                    if verbosity > 0:
                        mout.progress(total_price,budget,prepend='random BBS',append=f'#bbs={len(bb_set)} #prods={len(products)} ${total_price:.2f}')

                    old_bb_set = bb_set.copy()

                    iteration += 1

                if verbosity > 0:
                    mout.finish()

                # get products
                products = bb_set.get_products(all_compounds)
                mout.var('#products',len(products))
                for prod in products:
                    mout.out(f'{prod.name} {prod.smiles}')

                if not products:
                    # try again
                    mout.warning('No products, trying again')
                    return self.random_building_block_set(budget=budget, max_lead_time=max_lead_time, verbosity=verbosity, max_iter=max_iter)
                
                # remove bbs with no quantity
                count = bb_set.remove_unused(all_compounds)
                # mout.var('#bbs (removed)',count)

                mout.var('#bbs',len(bb_set))

                for bb in bb_set:

                    from rdkit.Chem import CanonSmiles
                    # if CanonSmiles(bb.smiles) != bb.smiles:
                    #     mout.error(f'Non-Canonical SMILES! {bb}')

                    if bb.required_amount:
                        mout.out(f'{bb.required_amount}mg of {bb}')
                    else:
                        mout.out(f'{mcol.error}{bb.required_amount}mg of {bb}')

                # plot_hist([bb.required_amount for bb in bb_set])
                
                total_price = bb_set.get_price()
                remaining = budget - total_price
                mout.var('total_price',f'${total_price:.2f}')
                mout.var('budget remaining',f'${remaining:.2f}')

                lead_times = bb_set.get_lead_times()
                this_max_lead_time = max(lead_times)
                mout.var('avg(lead_times)',f'{np.mean(lead_times):.1f}')
                mout.var('max(lead_times)',f'{this_max_lead_time:.0f}')

            elif by_set:

                assert isinstance(by_set,dict)
                assert not max_lead_time

                if verbosity:
                    mout.header(f'Choosing BBs up to ${budget} from {len(by_set)} product sets...')

                # shuffle the product sets
                for key in by_set:
                    random.shuffle(by_set[key])

                if start_with:
                    bb_set = start_with.copy()
                    products = bb_set.get_products(self.all_compounds)
                    products.immutable = False
                else:
                    bb_set = BuildingBlockSet()
                    products = CompoundSet(f'Products(BBS#{bbs_id})')

                iteration = -1
                old_bb_set = None

                while iteration < max_iter:

                    iteration += 1

                    if old_bb_set is not None:
                        start_price = old_bb_set.get_price()
                    else:
                        start_price = 0
                    
                    if verbosity:
                        mout.var("iteration",iteration)
                        mout.var("start_price",start_price)

                    if old_bb_set and iteration > 50:
                        if start_price > budget - 100:
                            if verbosity:
                                mout.error('Good enough')
                            break

                    # pick a random product set
                    try:
                        prod_set_key = random.choice(list(by_set.keys()))
                    except IndexError:
                        return None

                    # deal with empty set
                    if len(by_set[prod_set_key]) == 0:
                        del by_set[prod_set_key]
                        continue

                    # pick a random product
                    random_product = by_set[prod_set_key].pop()

                    # add the BBs
                    try:
                        for bb in random_product.building_blocks:
                            if not bb.purchaseable:
                                raise Exception(f'{bb} not purchaseable')
                            if bb not in bb_set:
                                bb_set.add(bb)
                                bb_set[bb].clear_amount()
                                if verbosity > 1:
                                    mout.out(f'Adding: {bb.smiles}')
                            bb_set[bb].increment_amount()
                    except Exception as e:
                        mout.error(str(e))
                        continue

                    # calculate total cost
                    total_price = bb_set.get_price()
                    remaining = budget - total_price
                    if verbosity > 1:
                        mout.var('total_price',total_price)
                        mout.var('budget remaining',remaining)

                    # lead times
                    if verbosity > 1:
                        lead_times = bb_set.get_lead_times()
                        this_max_lead_time = max(lead_times)
                        mout.var('avg(lead_times)',np.mean(lead_times))
                        if max_lead_time:
                            mout.var('max(lead_times)',this_max_lead_time,valCol=mcol.error if this_max_lead_time > max_lead_time else None)
                        else:
                            mout.var('max(lead_times)',this_max_lead_time)

                    if remaining < 0:
                        if verbosity > 1:
                            mout.error('Over budget')
                        if old_bb_set is None:
                            bb_set = BuildingBlockSet()
                            continue

                        bb_set = old_bb_set.copy()
                        continue

                    try:
                        products.add(random_product)
                    except ValueError:
                        continue

                    old_bb_set = bb_set.copy()

                
                # get products
                expected_products = bb_set.get_products(self.all_compounds)
                if debug:
                    mout.var('bb_set.total_bb_amount',bb_set.total_bb_amount)
                    mout.var('#products',len(products))
                    mout.var('#expected_products',len(expected_products))
                assert bb_set.total_bb_amount == 2*len(products)
                assert all([p.smiles in products.smiles for p in expected_products])

                bb_set._products = products
                bb_set._products.immutable = True

                if verbosity:

                    mout.var('#products',len(bb_set.products))
                    for prod in bb_set.products:
                        mout.out(f'{prod.name} {prod.smiles}')
                    
                    mout.var('#bbs',len(bb_set))

                    for bb in bb_set:
                        if bb.required_amount:
                            mout.out(f'{bb.required_amount}mg of {bb}')
                        else:
                            mout.out(f'{mcol.error}{bb.required_amount}mg of {bb}')
                
                    total_price = bb_set.get_price()
                    remaining = budget - total_price
                    mout.var('total_price',f'${total_price:.2f}')
                    mout.var('budget remaining',f'${remaining:.2f}')

                    lead_times = bb_set.get_lead_times()
                    this_max_lead_time = max(lead_times)
                    mout.var('avg(lead_times)',f'{np.mean(lead_times):.1f}')
                    mout.var('max(lead_times)',f'{this_max_lead_time:.0f}')

                bb_set.prepare_for_scoring(self.hit_features)

            else:

                available = CompoundSet('reminaining_products',self.all_products.compounds)
                available.shuffle()

                if debug:
                    print(f'{self.all_products=}')
                    print(f'{len(self.all_products)=}')

                    print(f'{available=}')
                    print(f'{len(available)=}')
                    
                    print(f'{len(self.all_compounds)=}')
                    print(f'{len(self.building_blocks)=}')

                if verbosity:
                    if max_lead_time:
                        mout.out(f'Choosing from {len(available)} product compounds up to ${budget} with lead_time <= {max_lead_time} days')
                        raise NotImplementedError
                    else:
                        mout.out(f'Choosing from {len(available)} product compounds up to ${budget}')

                bb_set = BuildingBlockSet()
                products = CompoundSet('from random BBS')

                iteration = -1
                old_bb_set = None
                over_budget_count = 0

                while iteration < max_iter:

                    iteration += 1

                    if old_bb_set is not None:
                        start_price = old_bb_set.get_price()
                    else:
                        start_price = 0
                    
                    if verbosity:
                        mout.var("iteration",iteration)
                        mout.var("start_price",start_price)

                    if old_bb_set and over_budget_count > retry_over_budget:
                        if start_price > budget - 100:
                            mout.header('Good enough')
                            break
                    
                    if len(available) < 1:
                        if verbosity:
                            mout.header('No more available')
                        break

                    # get a random product from available compounds
                    # random_product = random.choice(available)
                    # while random_product in products:
                        # random_product = random.choice(available)
                    random_product = available.pop()
                    if verbosity > 1:
                        mout.out(f'Trying: {random_product.name}')

                    # add the BBs
                    # try:

                    # skip if BBs too expensive
                    if bb_price_cutoff and any([bb.get_price(1) > bb_price_cutoff for bb in random_product.building_blocks]):
                        continue

                    for bb in random_product.building_blocks:
                        if not bb.purchaseable:
                            raise Exception(f'{bb} not purchaseable')
                        
                        if bb not in bb_set:
                            bb_set.add(bb)
                            bb_set[bb].clear_amount()
                            if verbosity > 1:
                                mout.out(f'Adding: {bb.smiles}')

                        if verbosity > 1:
                            mout.out(f'Incrementing: {bb.smiles}')
                        bb_set[bb].increment_amount()
                    
                    if verbosity > 1:
                        mout.var('bb_set.total_bb_amount',bb_set.total_bb_amount)
                        mout.var('#products-1',len(products))

                    # except Exception as e:
                        # mout.error(str(e))
                        # continue

                    # calculate total cost
                    total_price = bb_set.get_price()
                    remaining = budget - total_price
                    if verbosity > 1:
                        mout.var('total_price',total_price)
                        mout.var('budget remaining',remaining)

                    # lead times
                    if verbosity > 1:
                        lead_times = bb_set.get_lead_times()
                        this_max_lead_time = max(lead_times)
                        mout.var('avg(lead_times)',np.mean(lead_times))
                        if max_lead_time:
                            mout.var('max(lead_times)',this_max_lead_time,valCol=mcol.error if this_max_lead_time > max_lead_time else None)
                        else:
                            mout.var('max(lead_times)',this_max_lead_time)

                    if remaining < 0:
                        over_budget_count += 1
                        if verbosity > 1:
                            mout.error('Over budget')
                        if old_bb_set is None:
                            bb_set = BuildingBlockSet()
                        else:
                            bb_set = old_bb_set.copy()
                        continue

                    old_bb_set = bb_set.copy()

                    if verbosity > 1:
                        mout.out(f'Adding: {random_product.name}')
                    products.add(random_product)

                else:
                    mout.warning(f'Reached {max_iter=}')

                # get products
                # expected_products = bb_set.get_products(self.all_products)
                if debug:
                    mout.var('bb_set.total_bb_amount',bb_set.total_bb_amount)
                    mout.var('#products',len(products))
                    # mout.var('#expected_products',len(expected_products))
                assert bb_set.total_bb_amount == 2*len(products)
                # assert all([p.smiles in products.smiles for p in expected_products])
                # assert len(products) == len(expected_products)

                # # get products
                # products = bb_set.get_products(self.all_compounds)

                products.immutable = True
                bb_set._products = products

                if verbosity:

                    mout.var('#products',len(products))
                    for prod in products:
                        mout.out(f'{prod.name} {prod.smiles}')
                    
                    mout.var('#bbs',len(bb_set))

                    for bb in bb_set:
                        if bb.required_amount:
                            mout.out(f'{bb.required_amount}mg of {bb}')
                        else:
                            mout.out(f'{mcol.error}{bb.required_amount}mg of {bb}')
                
                    total_price = bb_set.get_price()
                    remaining = budget - total_price
                    mout.var('total_price',f'${total_price:.2f}')
                    mout.var('budget remaining',f'${remaining:.2f}')

                    lead_times = bb_set.get_lead_times()
                    this_max_lead_time = max(lead_times)
                    mout.var('avg(lead_times)',f'{np.mean(lead_times):.1f}')
                    mout.var('max(lead_times)',f'{this_max_lead_time:.0f}')

                bb_set.prepare_for_scoring(self.hit_features)

        return bb_set

    def fit_umap(self, n_dims=2, plot=True):

        if self.verbosity > 1:
            mout.debug('HIPPO.fit_umap()')

        import umap

        self._umap_df = self.compound_df

        self._umap_dims = n_dims

        self._reducer = umap.UMAP(n_components=self._umap_dims)

        # ignore metadata
        fingerprint_data = self._umap_df.loc[:, ~self._umap_df.columns.isin(self._comp_df_meta_keys)].values

        # fit a normaliser to the input data
        from sklearn.preprocessing import StandardScaler
        self._scaler = StandardScaler().fit(fingerprint_data)

        # normalise the input data
        scaled_data = self._scaler.transform(fingerprint_data)

        # fit the UMAP reducer to the protein data
        mout.debug('fitting scaled_data')
        self._transformer = self._reducer.fit(scaled_data)

        # get the transformed coordinates
        self._embedding = self._transformer.embedding_

        self._umap_df['u'] = self._embedding[:,0]
        if self._umap_dims == 2:
            self._umap_df['v'] = self._embedding[:,1]

        if plot:
            self.plot_umap()

    def plot_umap(self):

        assert self._umap_df is not None

        ### PLOTTING
        if self._umap_dims == 1:

            fig = px.scatter(
                self._umap_df,
                x='u',
                y='set_name',
                color='set_name',
                # symbol='set_name',
                hover_data=['name','smiles','set_name','is_pains','lead_time','cost_range_str','u']
            )

        elif self._umap_dims == 2:

            fig = go.Figure()

            px_fig = px.scatter(
                self._umap_df,
                x='u',
                y='v',
                color='set_name',
                hover_data=['name','smiles','set_name','is_pains','lead_time','cost_range_str','u','v']
            )

            for trace in px_fig.data:
                fig.add_trace(trace)

        else:
            raise Exception("multiple dimensions not implemented")

        # write the figure
        mp.write(f'umap_{self._umap_dims}d.html',fig)

    def export_building_blocks(self,path,cols):

        assert self.has_all_compounds

        all_bb_dicts = []

        for bb in self.building_blocks:

            bb_dict = { key:getattr(bb, key) for key in cols }

            all_bb_dicts.append(bb_dict)

        bb_df = pd.DataFrame(all_bb_dicts)
        bb_df.to_csv(open(path,'wt'))

    def cull_non_purchaseable_building_blocks(self, max_lead_time=None):

        assert self.has_purchase_info

        if self.has_non_purchaseable_been_culled:
            mout.warning(f'BBs have already been culled')

        self._building_blocks = BuildingBlockSet([bb for bb in self.building_blocks if bb.purchaseable])

        if self.verbosity:
            mout.var(f'#building_blocks (purchaseable)',f'{len(self.building_blocks)}')

        if max_lead_time is not None:

            self._building_blocks = BuildingBlockSet([bb for bb in self.building_blocks if bb.min_lead_time is not None and bb.min_lead_time < max_lead_time])
        
            if self.verbosity:
                mout.var(f'#building_blocks (min_lead_time: OK)',f'{len(self.building_blocks)}')

        mout.var('cost of all BBs (10mg)',self.building_blocks.get_price(10))
        
        lead_times = self.building_blocks.get_lead_times(10)

        mout.var('max(lead_times) (10mg)',max(lead_times))
        mout.var('avg(lead_times) (10mg)',np.mean(lead_times))

        self._has_non_purchaseable_been_culled = True

    def determine_all_products(self, cull_non_purchaseable):

        assert self.has_non_purchaseable_been_culled

        if self.has_products:
            mout.error(f'Products already determined')

        self.all_products = self.building_blocks.get_products(self.all_compounds)

        if cull_non_purchaseable:
            count = 0
            overlap = 0
            for comp_set in self.product_sets:
                for compound in comp_set:
                    if compound not in self.all_products:
                        comp_set.remove(compound)
                        count += 1

            if self.verbosity:
                mout.var('#non-product compounds removed',count)

        if self.verbosity:
            mout.var('#possible products',self.all_products.num_compounds)

        self._has_products = True

    def fingerprint_hits(self):

        assert self.has_hit_directory

        if self.has_hit_fingerprints:
            mout.error(f'Hits already fingerprinted')

        self._protein_feature_strs = self.generate_fingerprints(self.compound_sets[0])

        self._has_hit_fingerprints = True

    def fingerprint_products(self):

        assert self.has_products

        if self.has_product_fingerprints:
            mout.error(f'Products already fingerprinted')

        self._protein_feature_strs = self.generate_fingerprints(self.all_products)

        self._has_product_fingerprints = True

    def add_purchase_info(
        self,
        csv,
        smiles_col='Query SMILES',
        quoted_id_col='Quoted Mcule ID',
        lead_time_col='lead-time-days',
        amount_col='Quoted Amount (mg)',
        amount_unit='mg',
        price_col='Product price (USD)',
        price_unit='USD',
        lead_time_unit='days'
    ):
        
        if self.verbosity:
            mout.out(f'Getting purchase info from: {mcol.file}{csv}')

        purchase_df = pd.read_csv(csv)

        missing = []

        for index,row_df in purchase_df.iterrows():

            smiles = str(row_df[smiles_col])
            lead_time = float(row_df[lead_time_col])
            amount = float(row_df[amount_col])
            price = float(row_df[price_col])

            if smiles == 'nan':
                continue

            try:
                assert not np.isnan(lead_time), (csv, smiles)
                assert not np.isnan(amount), (csv, smiles)
                assert not np.isnan(price), (csv, smiles)
            except AssertionError as e:
                if self.verbosity > 2:
                    mout.warning(f'Skipping purchase row containing NaN: {e}')
                continue

            bb = self.building_blocks[smiles]

            if bb is None:
                missing.append(smiles)
                continue

            bb.add_purchase_info(
                lead_time=lead_time,
                price=price,
                amount=amount,
                row_df=row_df,
                amount_unit=amount_unit,
                price_unit=price_unit,
                lead_time_unit=lead_time_unit,
            )

        if missing:
            mout.warning(f'{len(missing)} unmatched BBs in purchase info')

    def generate_purchase_interpolators(self):
        
        for bb in list(self.building_blocks):

            if not bb.purchaseable:
                continue

            if not bb.has_purchase_interpolators:
                bb.generate_purchase_interpolators()

    def create_bb_map(self):

        assert self.has_all_compounds

        self._bb_id_to_smiles = {}
        self._bb_smiles_to_id = {}

        for i,bb in enumerate(self.building_blocks):
            bb.id = i
            self._bb_id_to_smiles[i] = bb.smiles
            self._bb_smiles_to_id[bb.smiles] = i

    def bb_id_to_smiles(self,bb_id):
        assert self.has_bb_map
        return self._bb_id_to_smiles[bb_id]

    def bb_smiles_to_id(self,smiles):
        assert self.has_bb_map
        return self._bb_smiles_to_id[smiles]

    def get_bb_set_hex_id(self,bb_set):
        fingerprint = []
        for i,s in self._bb_id_to_smiles.items():
            if s in bb_set:
                fingerprint.append(True)
            else:
                fingerprint.append(False)

        bin_str = '0b' + ''.join(['1' if x else '0' for x in fingerprint])
        number = eval(bin_str)
        return hex(number)

    def summary(self):
        mout.header(f'{self}')        

        mout.var('#building_blocks',len(self.building_blocks))
        mout.var('#compound_sets',len(self.compound_sets))
        mout.var('#compounds',self.num_compounds)

        mout.underline('compound sets:')
        for comp_set in self.compound_sets:
            mout.out(comp_set)

def plot_hist(x,show=True,bin_size=1):
    fig = go.Figure()
    trace = go.Histogram(x=x,xbins=dict(size=bin_size))
    fig.add_trace(trace)
    if show:
        fig.show()
    return fig
