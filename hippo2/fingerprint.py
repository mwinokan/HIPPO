
import mout
from pprint import pprint
import numpy as np

from molparse.rdkit.features import FEATURE_FAMILIES, COMPLEMENTARY_FEATURES
import molparse as mp

from collections import UserDict

CUTOFF_PADDING = 1.0

FEATURE_PAIR_CUTOFFS = {
    'Donor Acceptor': 3.5 + CUTOFF_PADDING,
    'Acceptor Donor': 3.5 + CUTOFF_PADDING,
    'NegIonizable PosIonizable': 4.5 + CUTOFF_PADDING,
    'PosIonizable NegIonizable': 4.5 + CUTOFF_PADDING,
    'Aromatic PosIonizable': 4.5 + CUTOFF_PADDING,
    'PosIonizable Aromatic': 4.5 + CUTOFF_PADDING,
    'Aromatic Aromatic': 6.0 + CUTOFF_PADDING,
    'Hydrophobe Hydrophobe': 4.5 + CUTOFF_PADDING,
}

FEATURE_METADATA = {}

class Fingerprint(UserDict):

    ### DUNDERS

    def __init__(self, pose, ligand_mol, protein_system):

        self.pose = pose
        self.data = self.calculate(ligand_mol, protein_system)

    def __repr__(self):
        return f'Fingerprint(#features={len(self)}, #interactions={self.num_interactions})'

    ### PROPERTIES

    @property
    def num_features(self):
        return len(self)

    @property
    def num_interactions(self):
        return sum([n for n in self.data.values()])
    
    ### METHODS

    def trim_keys(self, shared_keys):

        new_data = {}
        for key in shared_keys:
            new_data[key] = self.data[key]

        self.data = new_data

    def calculate(self, ligand_mol, protein_system):

        fingerprint = {}

        protein_features = protein_system.get_protein_features()
        comp_features = mp.rdkit.features_from_mol(ligand_mol)

        comp_features_by_family = {}
        for family in FEATURE_FAMILIES:
            comp_features_by_family[family] = [f for f in comp_features if f.family == family]

        for prot_feature in protein_features:

            prot_family = prot_feature.family

            complementary_family = COMPLEMENTARY_FEATURES[prot_family]

            complementary_comp_features = comp_features_by_family[complementary_family]

            cutoff = FEATURE_PAIR_CUTOFFS[f'{prot_family} {complementary_family}']

            valid_features = [f for f in complementary_comp_features if np.linalg.norm(prot_feature - f) <= cutoff]

            feature_key = prot_feature.family_name_number_chain_atoms_str

            if feature_key not in FEATURE_METADATA:

                FEATURE_METADATA[feature_key] = dict(
                    family=prot_feature.family,
                    atom_numbers=prot_feature.atom_numbers,
                    res_name=prot_feature.res_name,
                    res_number=prot_feature.res_number,
                    res_chain=prot_feature.res_chain,
                )

            fingerprint[feature_key] = len(valid_features)

        return fingerprint