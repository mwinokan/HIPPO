
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

        self._protein_feature_strs = None

    ### FACTORIES

    @classmethod
    def from_pickle(self, path):
        mout.debug('HIPPO.from_pickle()')
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

        tags = tags or []

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
        )

        self._compound_sets.append(comp_set)
        self._hit_compounds = self._compound_sets[-1]

        mout.success(f'Loaded {comp_set.num_compounds} compounds as "{comp_set.name}" ({self.hits.num_poses} poses)')

    def write_pickle(self,file):
        import molparse as mp
        mp.write(file,self)

    ### QUERIES

    def get_poses(self, tag, search=False):
        if search:
            return [p for p in self.all_poses if any([tag in t for t in p.tags])]
        else:
            return [p for p in self.all_poses if tag in p.tags]

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

    ### PLOTTING

    def plot_interaction_histogram(self, poses, subtitle=None):

        import plotly.express as px
        from .fingerprint import FEATURE_METADATA

        df = self._fingerprint_df(poses)

        plot_data = []
        
        for key in df.columns:
            count = int(df[key].sum())

            if not count:
                continue
            
            data = dict(str=key,count=count)

            data['family'] = FEATURE_METADATA[key]['family']
            data['res_name'] = FEATURE_METADATA[key]['res_name']
            data['res_number'] = FEATURE_METADATA[key]['res_number']
            data['res_chain'] = FEATURE_METADATA[key]['res_chain']
            data['atom_numbers'] = FEATURE_METADATA[key]['atom_numbers']

            data['res_name_number_chain_str'] = f"{FEATURE_METADATA[key]['res_name']} {FEATURE_METADATA[key]['res_number']} {FEATURE_METADATA[key]['res_chain']}"
        
            plot_data.append(data)

        plot_df = pd.DataFrame(plot_data)
        plot_df.sort_values(['res_chain', 'res_number', 'family'], inplace=True)
        plot_df

        fig = px.bar(plot_df, x='res_name_number_chain_str', y='count', color='family', hover_data=plot_df.columns)

        title='Leveraged protein features'

        if subtitle:
            fig.update_layout(title=f'{title} ({subtitle})')

        fig.update_layout(xaxis_title='Residue')
        fig.update_layout(yaxis_title='#Interactions')

        return fig

    def plot_interaction_punchcard(self, poses, subtitle, opacity=1.0):

        import plotly
        import plotly.express as px
        from .fingerprint import FEATURE_METADATA

        plot_data = []

        for pose in poses:

            fingerprint = pose._fingerprint

            # loop over each interaction in the pose
            for key, value in fingerprint.items():
                if not value:
                    continue
                
                data = dict(str=key,count=value)

                data['family'] = FEATURE_METADATA[key]['family']
                data['res_name'] = FEATURE_METADATA[key]['res_name']
                data['res_number'] = FEATURE_METADATA[key]['res_number']
                data['res_chain'] = FEATURE_METADATA[key]['res_chain']
                data['atom_numbers'] = FEATURE_METADATA[key]['atom_numbers'] 
                data['pose'] = pose.longname 

                data['res_name_number_chain_str'] = f"{FEATURE_METADATA[key]['res_name']} {FEATURE_METADATA[key]['res_number']} {FEATURE_METADATA[key]['res_chain']}"

                plot_data.append(data)

        plot_df = pd.DataFrame(plot_data)

        fig = px.scatter(plot_df, x='res_name_number_chain_str', y='family',marginal_x='histogram',marginal_y='histogram', hover_data=plot_df.columns, color='pose', title='Interaction Punch-Card')

        for trace in fig.data:
            if type(trace) == plotly.graph_objs._histogram.Histogram:
                trace.opacity = 1
                trace.xbins.size = 1
            else:
                trace['marker']['size'] = 10
                trace['marker']['opacity'] = opacity

        fig.update_layout(barmode='stack')
        fig.update_layout(scattermode='group', scattergap=0.75)

        return fig

    ### INTERNAL METHODS

    def _generate_fingerprints(self, poses):

        n = len(poses)

        fail_count = 0
        for i,pose in enumerate(poses):

            mout.progress(i, n, prepend='fingerprints', fill='#')

            try:
                pose.fingerprint

            except FailedToAssignBondOrders:
                mout.error(f'Failed to assign bond orders for {pose}')
                pose._fingerprint = None
                fail_count += 1
            except Exception as e:
                mout.error(e)
                pose._fingerprint = None
                fail_count += 1

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
    def all_tags(self):

        all_tags = []

        for thing in self.all_compounds.compounds + self.all_poses:
            for tag in thing.tags:
                if tag not in all_tags:
                    all_tags.append(tag)

        return all_tags
    
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
