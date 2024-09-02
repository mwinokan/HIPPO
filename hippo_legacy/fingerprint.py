import mout
from pprint import pprint
import numpy as np

from molparse.rdkit.features import FEATURE_FAMILIES, COMPLEMENTARY_FEATURES
import molparse as mp

from collections import UserDict

FEATURE_PAIR_CUTOFFS = {
    "Donor Acceptor": 3.5,
    "Acceptor Donor": 3.5,
    "NegIonizable PosIonizable": 4.5,
    "PosIonizable NegIonizable": 4.5,
    "Aromatic PosIonizable": 4.5,
    "PosIonizable Aromatic": 4.5,
    "Aromatic Aromatic": 6.0,
    "Donor Acceptor": 6.0,
    "Acceptor Donor": 6.0,
    "Hydrophobe Hydrophobe": 4.5,
}


class Fingerprint(UserDict):

    ### DUNDERS

    def __init__(self, compound):

        self.data = self.calculate(compound)

    ### PROPERTIES

    ### METHODS

    def trim_keys(self, shared_keys):

        new_data = {}
        for key in shared_keys:
            new_data[key] = self.data[key]

        self.data = new_data

    def calculate(self, compound):

        fingerprint = {}

        if compound._bound_pdb is not None:

            protein_features = compound.bound_system.get_protein_features()
            comp_features = compound.ligand_features

        else:

            protein_features = compound.protein_system.get_protein_features()
            comp_features = mp.rdkit.features_from_mol(compound.mol)

        comp_features_by_family = {}
        for family in FEATURE_FAMILIES:
            comp_features_by_family[family] = [
                f for f in comp_features if f.family == family
            ]

        for prot_feature in protein_features:

            prot_family = prot_feature.family

            complementary_family = COMPLEMENTARY_FEATURES[prot_family]

            complementary_comp_features = comp_features_by_family[complementary_family]

            cutoff = FEATURE_PAIR_CUTOFFS[f"{prot_family} {complementary_family}"]

            valid_features = [
                f
                for f in complementary_comp_features
                if np.linalg.norm(prot_feature - f) <= cutoff
            ]

            fingerprint[prot_feature.family_name_number_chain_atoms_str] = len(
                valid_features
            )
            # fingerprint[prot_feature.family_name_number_chain_str] = len(valid_features)

        return fingerprint

        # group = mp.rdkit.mol_to_AtomGroup(compound.mol)

        # group.plot3d(features=features)

        # exit()

        # raise Exception('Fingerprinting needs a bound PDB!')
