
import mout
import mcol

from .cset import CompoundSet
from .csetlist import CompoundSetList

from pathlib import Path
import pandas as pd
from pprint import pprint

class HIPPO:
    """Top-level HIPPO object

    Units
    -----

    * lead time: working days
    * reagant/product quantities: mg
    * pricing: arbitrary / whatever the quoting is

    """

    def __init__(self, name):

        self._name = name

        self._max_bb_price = None
        self._min_bb_quantity = None
        self._max_lead_time = None
        
        self._hit_compounds = None

        self._compound_sets = CompoundSetList()

    ### FACTORIES

    ### PUBLIC METHODS

    def set_cutoffs(self, max_lead_time, max_bb_price, min_bb_quantity):
        self._max_lead_time = max_lead_time
        self._max_bb_price = max_bb_price # do we need this?
        self._min_bb_quantity = min_bb_quantity

    def add_hits(self, metadata_path, pdb_path, pdb_pattern='**/*-x????_??_bound.pdb', tags=None):
        
        mout.debug(f'{self}.add_hits()')
        
        ### checks
        
        if 'hits' in [s.name for s in self.compound_sets]:
            mout.error(f'CompoundSet "hits" already exists')

        if not isinstance(metadata_path, Path):
            metadata_path = Path(metadata_path)

        if not isinstance(pdb_path, Path):
            pdb_path = Path(pdb_path)

        tags = tags or []

        ### metadata

        mout.var('metadata_path',str(metadata_path),valCol=mcol.file)

        try:
            metadata_df = pd.read_csv(metadata_path)
        except FileNotFoundError as e:
            mout.error(f'FileNotFoundError: {e}')
            return

        mout.header('metadata columns:')
        pprint(list(metadata_df.columns))

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
        )

        self._compound_sets.append(comp_set)
        self._hit_compounds = self._compound_sets[-1]

        mout.success(f'Loaded {comp_set.num_compounds} compounds as "{comp_set.name}" ({self.hits.num_poses} poses)')

    def summary(self):
        mout.header(f'{self}')  

        mout.var('max_lead_time', self.max_lead_time, unit='workdays')      
        mout.var('max_bb_price', self.max_bb_price, unit='$')      
        mout.var('min_bb_quantity', self.min_bb_quantity, unit='mg')      

        # mout.var('#building_blocks',len(self.building_blocks))
        mout.var('#compound_sets',len(self.compound_sets))
        mout.var('#compounds',self.num_compounds)
        mout.var('#poses',self.num_poses)

        if self.compound_sets:
            mout.out('')
            mout.underline('compound sets:')
            for comp_set in self.compound_sets:
                mout.out(comp_set)

    ### PRIVATE METHODS

    ### PROPERTIES

    @property
    def name(self):
        return self._name

    @property
    def num_compounds(self):
        return len(self.all_compounds)

    @property
    def all_compounds(self):
        if len(self.compound_sets) == 1:
            return self.compound_sets[0]
        return set().union(*self.compound_sets)

    @property
    def num_poses(self):
        return len(self.all_poses)

    @property
    def all_poses(self):
        return sum([c.poses for c in self.all_compounds], [])

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
    
    ### DUNDERS

    def __repr__(self):
        return f'HIPPO({self.name})'
