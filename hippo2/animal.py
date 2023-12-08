
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

from rdkit import Chem

from pprint import pprint

from .tools import df_row_to_dict

class HIPPO:

    """Top-level HIPPO object

    Units
    -----

    * lead time: working days
    * reagant/product quantities: mg
    * pricing: arbitrary / whatever the quoting is

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

    def __init__(self, project_name, target_name):

        self._name = project_name
        self._target_name = target_name

        self._max_bb_price = None
        self._min_bb_quantity = None
        self._max_lead_time = None
        
        self._hit_compounds = None

        self._compound_sets = CompoundSetList()
        self._building_blocks = CompoundSet('building_blocks')

        self._protein_feature_strs = None

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
        )

        if 'hits' in self.compound_sets:
            self._compound_sets['hits'] = comp_set
        else:
            self._compound_sets.append(comp_set)
        
        self._hit_compounds = self._compound_sets['hits']

        mout.success(f'Loaded {comp_set.num_compounds} compounds as "{comp_set.name}" ({self.hits.num_poses} poses)')

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

        # return meta_dict

    def add_elabs(self, mol_paths, pdb_paths, metadata, tags=None, set_name='elabs', append=False, overwrite=False, first_is_base=True):

        ### checks
        
        if not overwrite and set_name in self.compound_sets:
            mout.error(f'CompoundSet "{set_name}" already exists')
            return

        tags = tags or [set_name]
        tags = TagSet(tags)

        if append:
            comp_set = self.compound_sets[set_name]
        else:
            comp_set = CompoundSet(set_name)

        meta_df = pd.read_csv(metadata)

        for mol, pdb in zip(mol_paths, pdb_paths):

            if not isinstance(mol, Path):
                mol = Path(mol)

            comp_name = mol.name.replace('.minimised.mol','')

            compound = Compound.from_mol(comp_name, mol, tags=tags)

            pose = Pose.from_bound_pdb(compound, pdb)

            compound.add_pose(pose)

            comp_set.add(compound)

            ### inspirations

            ### base

            ### synthetic routes

        if not append:
            if overwrite and set_name in self.compound_sets:
                self.compound_sets[set_name] = comp_set
            else:
                self.compound_sets.append(comp_set)
            # self._base_compounds = self._compound_sets[-1]

        mout.success(f'Loaded {comp_set.num_compounds} compounds as "{comp_set.name}" ({comp_set.num_poses} poses)')

        return meta_df

    def summary(self):
        mout.header(f'{self}')  

        mout.var('target_name', self.target_name)      
        mout.var('max_lead_time', self.max_lead_time, unit='workdays')      
        mout.var('max_bb_price', self.max_bb_price, unit='$')      
        mout.var('min_bb_quantity', self.min_bb_quantity, unit='mg')      

        # mout.var('#building_blocks',len(self.building_blocks))
        mout.var('#compound_sets',len(self.compound_sets))
        mout.var('#compounds',self.num_compounds)
        mout.var('#poses',self.num_poses)
        if self.num_poses:
            mout.var('#tags',self.num_tags)

        if self.compound_sets:
            mout.out('')
            mout.underline('compound sets:')
            for comp_set in self.compound_sets:
                mout.out(comp_set)
        
            all_tags = self.all_tags
            if all_tags:
                mout.out('')
                mout.underline('tags:')
                for tag in all_tags:

                    num_compounds = len(self.get_compounds(tag))
                    num_poses = len(self.get_poses(tag))

                    mout.out(f'{tag} #compounds={num_compounds}, #poses={num_poses}')

    def write_pickle(self,file):
        import molparse as mp
        mp.write(file,self)

    ### QUERIES

    def get_compounds(self, tag, search=False):
        if search:
            return CompoundSet(f'{tag} in tags (search)',[c for c in self.all_compounds if any([tag in t for t in c.tags])])
        else:
            return CompoundSet(f'{tag} in tags',[c for c in self.all_compounds if tag in c.tags])

    def get_poses(self, tag, search=False):
        if search:
            return PoseSet([p for p in self.all_poses if any([tag in t for t in p.tags])])
        else:
            return PoseSet([p for p in self.all_poses if tag in p.tags])

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

    def plot_tag_statistics(self):
        from .plotting import plot_tag_statistics
        return plot_tag_statistics(self)

    def plot_interaction_histogram(self, poses, subtitle=None):
        from .fingerprint import FEATURE_METADATA
        from .plotting import plot_interaction_histogram
        return plot_interaction_histogram(self, poses, FEATURE_METADATA, subtitle)

    def plot_interaction_punchcard(self, poses, subtitle, opacity=1.0):
        from .fingerprint import FEATURE_METADATA
        from .plotting import plot_interaction_punchcard
        return plot_interaction_punchcard(self, poses, FEATURE_METADATA, subtitle, opacity)

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
    def _add_compound_metadata(self, compound, meta_dict):

        ### inspirations

        if 'hit_names' in meta_dict:
            for hit_name in meta_dict['hit_names'].split():
                hit_name = hit_name.replace(f'{self.target_name}-','')
                compound.inspirations.append(self.get_hit_pose(hit_name))

        else:
            mout.error(f'No inspiration data for {compound}')

        ### synthetic routes

        num_reactants = len([k for k in meta_dict if 'reactant' in k and 'smi' in k]) ### OLD FORMAT

        # pprint(meta_dict)

        # mout.debug(num_reactants)

        if num_reactants == 0:
            return

        elif num_reactants == 2:

            # mout.debug('Adding single-step reaction')

            # grab metadata
            mout.warning('Reaction type is unknown')
            reaction_type = 'Unknown'

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

            # create the reaction
            reaction = Reaction(reaction_type, compound, [reactant1, reactant2], product_amount=amounts, amounts=amounts)

            compound.add_reaction(reaction)

        elif num_reactants > 2:
            mout.error('Not yet implemented')
            raise NotImplementedError

        else:
            mout.error(f'reactants for {compound} = {num_reactants}')
            raise Exception('Unsupported number of reactants')

    def _get_or_create_building_block(self, smiles):

        # N.B. BuildingBlock.smiles == BuildingBlock.name

        if smiles in self.building_blocks:
            return self.building_blocks[smiles]

        else:
            bb = BuildingBlock(smiles)
            self.building_blocks.add(bb)
            return bb

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
        return PoseSet(sum([c.poses for c in self.all_compounds], PoseSet()))

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
        return self.compound_sets['bases']
    
    @property
    def target_name(self):
        return self._target_name

    @property
    def building_blocks(self):
        return self._building_blocks
    
    ### DUNDERS

    def __repr__(self):
        return f'HIPPO({self.name})'

class FailedToAssignBondOrders(Exception):
    pass
