
import mout
import mcol
import molparse as mp

import pandas as pd

from .set import CompoundSet
from .compound import Compound
from .block import BuildingBlock

import pickle

from pprint import pprint

from rdkit.Chem import PandasTools

import plotly.graph_objects as go

from .compound import FailedToAssignBondOrders

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

    def __init__(self, project_key, verbosity=3):
        self._project_key = project_key
        self._verbosity = verbosity

        self._protein_reference_path = None
        self._protein_system = None
        self._hit_pdbs = None
        self._metadata_df = None
        self._compound_sets = []
        self._missed_features = None
        self._building_blocks = {}

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
        return sum([s.compounds for s in self.compound_sets],[])

    @property
    def fingerprinted_compounds(self):
        return [c for c in self.all_compounds if c.fingerprint is not None]

    @property
    def num_compounds(self):
        return len(self.all_compounds)

    @property
    def num_building_blocks(self):
        return len(self.building_blocks)

    @property
    def building_blocks(self):
        return self._building_blocks

    ### METHODS

    def write_pickle(self,file):
        mp.write(file,self)
    
    def add_protein_reference(self, path):

        if self.verbosity > 1:
            mout.debug('HIPPO.add_protein_reference()')
            mout.var('path',str(path),valCol=mcol.file)

        self._protein_reference_path = path
        self._protein_system = mp.parsePDB(path,verbosity=self.verbosity-2)

        self.protein_system.name = f'{self.project_key} [protein reference]'

        if self.verbosity > 1:
            self.protein_system.summary()

        for chain in self.protein_system.chains:
            assert chain.type == 'PRO'

        if self.verbosity:
            mout.success(f'Parsed protein PDB')

    def add_hit_metadata(self, path):

        if self.verbosity > 1:
            mout.debug('HIPPO.add_hit_metadata()')
            mout.var('path',str(path),valCol=mcol.file)

        self._metadata_df = pd.read_csv(path)

        if self.verbosity > 2:
            mout.header('metadata columns:')
            pprint(list(self.metadata_df.columns))

        if self.verbosity:
            mout.success(f'Parsed metadata CSV')

    def add_hit_directory(self, name='hits', path=None, pattern='*-x????_??_bound.pdb'):

        assert path is not None

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

        return 

    def add_compound_from_sdf(self, name, path, idName='name',molColName='mol'):

        if self.verbosity > 1:
            mout.debug('HIPPO.add_compound_sdf()')
            mout.var('name',name,valCol=mcol.arg)
            mout.var('path',str(path),valCol=mcol.file)
            mout.var('idName',idName,valCol=mcol.arg)
            mout.var('molColName',molColName,valCol=mcol.arg)
            
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

    def add_product_compounds(self, name, metadata, data_path, mol_pattern):

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

        for index, (smiles, mol_weight, metadata) in meta_df[['reactant1_smi','reactant1_mw','reactant1_metadata']].iterrows():
            
            if smiles not in self._building_blocks:
                self._building_blocks[smiles] = BuildingBlock(smiles, mol_weight, metadata)

        for index, (smiles, mol_weight, metadata) in meta_df[['reactant2_smi','reactant2_mw','reactant2_metadata']].iterrows():
            
            if smiles not in self._building_blocks:
                self._building_blocks[smiles] = BuildingBlock(smiles, mol_weight, metadata)

        if self.verbosity:
            mout.success(f'Found {self.num_building_blocks} building blocks for "{name}"')

        comp_set = CompoundSet(name=name)

        mols = list(data_path.glob(mol_pattern))

        progress = len(mols) > 100

        for i,mol_file in enumerate(mols):

            name = mol_file.name.split('.')[0]

            if progress:
                mout.progress(i,len(mols),prepend=f'Loading compounds "{name}"',append=name,fill='🦛')

            mol = mp.parse(mol_file, verbosity=0)

            compound = Compound.from_rdkit_mol(name, mol)

            if 'base' in name:
                place_type = 'base'
            elif 'frag':
                place_type = 'frag'
            else:
                raise Exception

            meta_row = meta_df[meta_df[f'{place_type}_placed_name'] == name]
            assert len(meta_row)
            
            minimisation_failed = meta_row[f'{place_type}_fail_minimization_without_constraint'].values[0]

            if minimisation_failed:
                mout.error(f'Skipping b/c minimisation failed: {compound}')
                continue

            compound._smiles = meta_row['product_smi'].values[0]
            compound._is_pains = meta_row['is_pains'].values[0]
            compound._protein_system = self.protein_system

            comp_set.add_compound(compound)
        
        mout.finish()

        self._compound_sets.append(comp_set)

        if self.verbosity:
            mout.success(f'Loaded {comp_set.num_compounds} compounds "{comp_set.name}"')

    def generate_fingerprints(self):

        if self.verbosity > 1:
            mout.debug('HIPPO.generate_fingerprints()')

        fail_count = 0
        for i,compound in enumerate(self.all_compounds):
            if self.verbosity:
                mout.progress(i, self.num_compounds, prepend='fingerprints', append=f'{mcol.bold}{compound.set_name}{mcol.clear} {compound.name}', fill='🦛')

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

        self._protein_feature_strs = set(compound.fingerprint.keys()).intersection(*[set(c.fingerprint.keys()) for c in self.fingerprinted_compounds])

        for compound in self.fingerprinted_compounds:
            compound.fingerprint.trim_keys(self._protein_feature_strs)

        if self.verbosity > 1:
            mout.var('#feature_dims',len(self._protein_feature_strs))

        if self.verbosity:
            mout.finish(append='')
            if fail_count:
                mout.error(f"Failed to fingerprint {fail_count} compounds")
            mout.success(f'Fingerprinted {self.num_compounds - fail_count}/{self.num_compounds} compounds')

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

    def compare_feature_sets(self, name1, name2, set1, set2, print_diff=False):

        in_1_but_not_2 = set1.difference(set2)
        mout.var(f'# only in {name1}',len(in_1_but_not_2))
        mout.var(f'% missed {name1} features',f'{len(in_1_but_not_2)/len(set1):.1%}')
        mout.var(f'% covered {name1} features',f'{(len(set1)-len(in_1_but_not_2))/len(set1):.1%}')
        if print_diff:
            pprint(in_1_but_not_2)

        in_2_but_not_1 = set2.difference(set1)  
        mout.var(f'# only in {name2}',len(in_2_but_not_1))
        mout.var(f'% missed {name1} features',f'{len(in_2_but_not_1)/len(set2):.1%}')
        mout.var(f'% covered {name1} features',f'{(len(set2)-len(in_2_but_not_1))/len(set2):.1%}')
        if print_diff:
            pprint(in_2_but_not_1)

        return in_1_but_not_2, in_2_but_not_1
