
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
        # self._max_lead_time = None

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

            data.append(comp.as_dict)

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

    def add_reaction_products(self, name, metadata, data_path, mol_pattern, skip_unconstrained_minimisation_fail=False, expecting=None):

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

        progress = len(mols) > 100

        for i,mol_file in enumerate(mols):

            prod_name = mol_file.name.split('.')[0]

            if progress:
                mout.progress(i,len(mols),prepend=name,append=prod_name,fill='#')

            mol = mp.parse(mol_file, verbosity=0)

            compound = Compound.from_rdkit_mol(prod_name, mol)

            if 'base' in prod_name:
                place_type = 'base'
            elif 'frag':
                place_type = 'frag'
            else:
                raise Exception

            if f'{place_type}_placed_name' in meta_df:
                meta_row = meta_df[meta_df[f'{place_type}_placed_name'] == prod_name]
                assert len(meta_row)
            else:
                meta_row = meta_df[meta_df['name'] == prod_name]
                assert len(meta_row)
            
            if f'{place_type}_fail_minimization_without_constraint' in meta_row:
                compound._minimization_failed_without_constraint = meta_row[f'{place_type}_fail_minimization_without_constraint'].values[0]

            if 'product_smi' in meta_row:
                compound._smiles = meta_row['product_smi'].values[0]
            elif 'smiles' in meta_row:
                compound._smiles = meta_row['smiles'].values[0]
            else:
                mout.error(f"Couldn't get SMILES from metadata {index=}")
                raise Exception("No SMILES")

            if 'is_pains' in meta_row:
                compound._smiles = meta_row['is_pains'].values[0]
            else:
                mout.error(f"Couldn't get PAINS from metadata {index=}")
                raise Exception("No PAINS")

            if '∆∆G' in meta_row:
                compound._fragmenstein_ddG = meta_row['∆∆G'].values[0]
            else:
                # get fragmenstein json
                try:
                    with open(str(mol_file).replace('.mol','.json')) as f:
                        mol_json_data = json.load(f)
                        # mout.json(mol_json_data)

                    compound._fragmenstein_ddG = mol_json_data['Energy']['ligand_ref2015']['total_score'] - mol_json_data['Energy']['unbound_ref2015']['total_score']
                    compound._fragmenstein_mRMSD = mol_json_data['mRMSD']
                    compound._fragmenstein_too_moved = compound._fragmenstein_mRMSD > 1.0
                    compound._fragmenstein_too_contorted = compound._fragmenstein_ddG > 0.0

                except FileNotFoundError:
                    pass

            if 'LE' in meta_row:
                compound._fragmenstein_ligand_efficiency = meta_row['LE'].values[0]

            if 'outcome' in meta_row:
                compound._fragmenstein_outcome = meta_row['outcome'].values[0]

            compound._protein_system = self.protein_system
            
            reactant1_smiles = meta_row['reactant1_smi'].values[0]
            reactant2_smiles = meta_row['reactant2_smi'].values[0]
            
            bb1 = self.building_blocks[reactant1_smiles]
            bb2 = self.building_blocks[reactant2_smiles]
            compound._building_blocks = BuildingBlockSet([bb1,bb2])

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

    def random_building_block_set(self, budget=None, size=None, max_lead_time=None, verbosity=1, max_iter=1000, start_with=None, top_down=False, debug=True):

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

            else:

                # get all produceable compounds
                available = self.all_products
                # mout.debug(available)
                # self.building_blocks.get_products(self.all_compounds)
                # available = self.all_products
                # mout.debug(available._elements)
                # all_bbs = self.building_blocks

                if debug:
                    print(f'{self.all_products=}')
                    print(f'{id(self.all_products)=}')
                    print(f'{len(self.all_products)=}')
                    print(f'{len(self.all_compounds)=}')
                    print(f'{len(self.building_blocks)=}')

                if max_lead_time:
                    mout.out(f'Choosing from {len(available)} product compounds up to ${budget} with lead_time <= {max_lead_time} days')
                else:
                    mout.out(f'Choosing from {len(available)} product compounds up to ${budget}')

                bb_set = BuildingBlockSet()
                products = CompoundSet('from random BBS')

                iteration = 0
                old_bb_set = None

                while iteration < max_iter:

                    # get a random product from available compounds
                    random_product = random.choice(available)
                    while random_product in products:
                        random_product = random.choice(available)
                    products.add(random_product)

                    # add the BBs
                    for bb in random_product.building_blocks:
                        if bb not in bb_set:
                            bb_set.add(bb)
                            bb_set[bb].required_amount = 0
                            if verbosity > 1:
                                mout.out(f'Adding: {bb.smiles}')
                        bb_set[bb].required_amount += 1

                    # calculate total cost
                    total_price = bb_set.get_price()
                    remaining = budget - total_price
                    if verbosity > 1:
                        mout.var('total_price',total_price)
                        mout.var('budget remaining',remaining)

                    # lead times
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
                        if old_bb_set is None:
                            # try again
                            return self.random_building_block_set(budget=budget, size=size, max_lead_time=max_lead_time, verbosity=verbosity, max_iter=max_iter, start_with=start_with, top_down=top_down)

                        bb_set = old_bb_set.copy()
                        break

                    # if verbosity > 0:
                    #     mout.progress(total_price,budget,prepend='random BBS',append=f'#bbs={len(bb_set)} #prods={len(products)} ${total_price:.2f}')

                    old_bb_set = bb_set.copy()

                    iteration += 1

                # get products
                products = bb_set.get_products(self.all_compounds)
                # import time
                # products.name = f'{products.name} {time.time()}'
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
            for comp_set in self.compound_sets[1:]:
                print(comp_set,len(comp_set))
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

            assert not np.isnan(lead_time), (csv, smiles, lead_time, amount, price)
            assert not np.isnan(amount), (csv, smiles, lead_time, amount, price)
            assert not np.isnan(price), (csv, smiles, lead_time, amount, price)

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
            mout.error(f'Could not identify {len(missing)} BBs by SMILES')

    def generate_purchase_interpolators(self):
        
        for bb in list(self.building_blocks)[1:]:

            if not bb.purchaseable:
                continue

            if not bb.has_purchase_interpolators:
                bb.generate_purchase_interpolators()

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
