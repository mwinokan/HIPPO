
import mout
import mcol
import molparse as mp

import pandas as pd

from .set import CompoundSet

from pprint import pprint

from rdkit.Chem import PandasTools

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
    def num_compounds(self):
        return len(self.all_compounds)

    ### METHODS
    
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

    def add_compound_sdf(self, name, path, idName='name',molColName='mol'):

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

        comp_set = CompoundSet.from_df(name=name, df=mol_df[[idName,molColName]])
        
        if self.verbosity > 3:
            print(comp_set)
            pprint(comp_set.compounds)

        self._compound_sets.append(comp_set)

        if self.verbosity:
            mout.success(f'Loaded {comp_set.num_compounds} compounds "{comp_set.name}"')

    def generate_fingerprints(self):

        if self.verbosity > 1:
            mout.debug('HIPPO.generate_fingerprints()')
        
        fail_count = 0
        for i,compound in enumerate(self.all_compounds):
            if self.verbosity:
                mout.progress(i, self.num_compounds, prepend='fingerprints', append=f'{mcol.bold}{compound.set_name}{mcol.clear} {compound.name}')

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


        if self.verbosity:
            mout.finish(append='')
            if fail_count:
                mout.error(f"Failed to fingerprint {fail_count} compounds")
            mout.success(f'Fingerprinted {self.num_compounds - fail_count}/{self.num_compounds} compounds')
