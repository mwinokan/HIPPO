
import os
from rdkit import RDConfig, Chem
from rdkit.Chem import AllChem
import mout
import molparse as mp
import numpy as np
from .fingerprint import Fingerprint

class Compound:

    _fdef = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef'))

    ### DUNDERS

    def __init__(self, name):

        self._name = name

        # all
        self._smiles = None
        self._fingerprint = None
        self._ligand_group = None
        self._set_name = None
        self._ligand_features = None
        self._is_pains = None
        self._building_blocks = None
        
        # from_rdkit_mol
        self._mol = None
        self._protein_system = None

        # from_bound_pdb
        self._bound_pdb = None
        self._bound_system = None
        self._pose_name = None
        self._crystal_name = None
        self._alternate_name = None
        self._site_name = None
        self._pdb_entry = None
        self._site_index = None
        self._chain_char = None

        # fragmenstein
        self._fragmenstein_ddG = None
        self._fragmenstein_ligand_efficiency = None
        self._fragmenstein_outcome = None

    def __repr__(self):
        if self._building_blocks:
            return f'Compound("{self.name}", "{self.smiles}", cost={self.cost_range_str}, lead_time={self.lead_time})'
        else:
            return f'Compound("{self.name}", "{self.smiles}")'

    ### FACTORIES

    @classmethod
    def from_bound_pdb(cls, name, pdb, metadata):

        self = cls.__new__(cls)

        self.__init__(name)

        self._chain_char = name[-1]
        self._site_index = int(name[-2])

        self._bound_pdb = pdb

        self._pose_name = metadata['crystal_name']
        self._crystal_name = metadata['RealCrystalName']
        self._smiles = metadata['new_smiles'] or metadata['smiles']
        self._alternate_name = metadata['alternate_name']
        self._site_name = metadata['site_name']
        self._pdb_entry = metadata['pdb_entry']

        return self
    
    @classmethod
    def from_rdkit_mol(cls, name, mol, protonate=False, verbosity=1):

        self = cls.__new__(cls)

        self.__init__(name)

        self._mol = mol

        if protonate:
            try:
                prot_mol = mp.rdkit.protonate(self.mol,verbosity=verbosity-1)
            except ValueError as e:
                mout.error(f'Cound not protonate: {e}')
                self._protonation_failed = True
            else:
                self._mol = prot_mol
                self._protonation_failed = False

        return self

    ### PROPERTIES

    @property
    def name(self):
        return self._name

    @property
    def set_name(self):
        return self._set_name

    @property
    def mol(self):
        if self._mol is None and self._bound_pdb:

            # get the coordinates from the bound PDB but use SMILES to solve bond orders, etc.

            # build the molecules
            m_smiles = Chem.MolFromSmiles(self.smiles)
            m_pdb = self.ligand_group.rdkit_mol

            # get bond orders from smiles
            try:
                m_pdb = AllChem.AssignBondOrdersFromTemplate(m_smiles, m_pdb)
            except ValueError as e:
                raise FailedToAssignBondOrders()

            # protonate
            m_pdb_prot = Chem.AddHs(m_pdb)

            # solve for a pose
            ps = AllChem.ETKDGv3()
            AllChem.EmbedMolecule(m_pdb_prot,ps)

            self._mol = m_pdb_prot

        return self._mol

    @property
    def is_pains(self):
        return self._is_pains

    @property
    def smiles(self):
        if self._smiles is None:
            self._smiles = mp.rdkit.mol_to_smiles(self.mol)
        return self._smiles

    @property
    def protein_system(self):
        return self._protein_system

    @property
    def fingerprint(self):
        return self._fingerprint

    @property
    def as_dict(self):

        data = {}

        data['name'] = self.name
        data['smiles'] = self.smiles
        data['set_name'] = self.set_name
        data['is_pains'] = self.is_pains
        data['cost_range_str'] = self.cost_range_str
        data['lead_time'] = self.lead_time
        
        for key, value in self.fingerprint.items():
            data[key] = value

        return data

    @property
    def ligand_group(self):
        if self._ligand_group is None:

            # if self._bound_pdb:

            lig_residues = self.bound_system['rLIG']
            lig_residues = [l for l in lig_residues if l.chain == self._chain_char]

            split_lig_residues = []
            for lig in lig_residues:
                split_lig_residues += lig.split_by_site()

            self._ligand_group = split_lig_residues[self._site_index]

            # else:

            #     try:
            #         self._mol = mp.rdkit.protonate(self.mol)
            #     except ValueError as e:
            #         mout.error(e)
            #         return None

            #     self._ligand_group = mp.rdkit.mol_to_AtomGroup(self.mol)

        return self._ligand_group

    @property
    def ligand_features(self):
        if self._ligand_features is None:

            assert self._bound_pdb

            # get the features
            raw_features = self._fdef.GetFeaturesForMol(self.mol)

            self._ligand_features = []
            for feat in raw_features:

                indices = feat.GetAtomIds()
                family = feat.GetFamily()

                # position from indices
                if len(indices) == 1:
                    position = self.ligand_group.atoms[indices[0]].np_pos
                else:
                    position = np.mean([self.ligand_group.atoms[i].np_pos for i in indices],axis=0)

                f_dict = dict(
                    family=family,
                    position=position,
                    indices=indices,
                    x=position[0],
                    y=position[1],
                    z=position[2],
                )

                feature = mp.rdkit.Feature(
                    family=f_dict['family'], 
                    atoms=[self.ligand_group.atoms[i] for i in indices],
                    position=position,
                    sidechain=None,
                    res_name='LIG',
                    res_number=None,
                    res_chain=self._chain_char,
                )

                self._ligand_features.append(feature)

        return self._ligand_features

    @property
    def bound_pdb(self):
        if self._bound_pdb is None:
            mout.error(f'{self} has no self._bound_pdb',fatal=True)
        return self._bound_pdb

    @property
    def bound_system(self):
        if self._bound_system is None:
            self._bound_system = mp.parsePDB(self.bound_pdb, verbosity=0, alternative_site_warnings=False)
        return self._bound_system

    @property
    def building_blocks(self):
        return self._building_blocks

    @property
    def cost_range_str(self):
        if self.building_blocks:
        
            cost_min = 0
            cost_max = 0
            for bb in self.building_blocks:

                if bb._cost_min is None:
                    return None
                if bb._cost_max is None:
                    return None

                # assert bb._cost_unit == '/g'
                cost_min += bb._cost_min
                cost_max += bb._cost_max

            return f'${cost_min}-{cost_max}'

        return None
    
    @property
    def lead_time(self):
        if self.building_blocks:
            lead_times = []
            for bb in self.building_blocks:
                if bb.lead_time is None:
                    return None
                lead_times.append(bb.lead_time)
            lead_time = max(lead_times)
            return lead_time
        return None
    
    ### METHODS

    def calculate_fingerprint(self):
        self._fingerprint = Fingerprint(self)

class FailedToAssignBondOrders(Exception):
    pass