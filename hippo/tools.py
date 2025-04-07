import re
import numpy as np
from molparse.rdkit import mol_from_smiles
from rdkit.Chem.inchi import MolToInchiKey
from rdkit.Chem import MolFromSmiles, MolToSmiles, AddHs, RemoveHs
import mcol
import mout
from datetime import datetime
from string import ascii_uppercase

import mrich


def df_row_to_dict(df_row):
    """

    :param df_row:

    """

    assert len(df_row) == 1, f"{len(df_row)=}"

    data = {}

    for col in df_row.columns:

        if col == "Unnamed: 0":
            continue

        value = df_row[col].values[0]

        if not isinstance(value, str) and np.isnan(value):
            value = None

        data[col] = value

    return data


def remove_other_ligands(sys, residue_number, chain):
    """

    :param sys:
    :param residue_number:
    :param chain:

    """

    ligand_residues = [r.number for r in sys["rLIG"] if r.number != residue_number]

    # if ligand_residues:
    for c in sys.chains:
        if c.name != chain:
            c.remove_residues(names=["LIG"], verbosity=0)
        elif ligand_residues:
            c.remove_residues(numbers=ligand_residues, verbosity=0)

    # print([r.name_number_str for r in sys['rLIG']])

    assert (
        len([r.name_number_str for r in sys["rLIG"]]) == 1
    ), f"{sys.name} {[r.name_number_str for r in sys['rLIG']]}"

    return sys


def inchikey_from_smiles(smiles):
    """

    :param smiles:

    """
    mol = mol_from_smiles(smiles)
    return MolToInchiKey(mol)


def flat_inchikey(smiles):
    """

    :param smiles:

    """
    smiles = sanitise_smiles(smiles)
    return inchikey_from_smiles(smiles)


def remove_isotopes_from_smiles(smiles):
    """

    :param smiles:

    """

    mol = MolFromSmiles(smiles)

    atom_data = [(atom, atom.GetIsotope()) for atom in mol.GetAtoms()]

    for atom, isotope in atom_data:
        if isotope:
            atom.SetIsotope(0)

    return MolToSmiles(mol)


def smiles_has_isotope(smiles, regex=True):
    """

    :param smiles:
    :param regex:  (Default value = True)

    """
    if regex:
        return re.search(r"([\[][0-9]+[A-Z]+\])", smiles)
    else:
        mol = MolFromSmiles(smiles)
        return any(atom.GetIsotope() for atom in mol.GetAtoms())


def sanitise_smiles(s, verbosity=False, sanitisation_failed="error", radical="error"):
    """

    :param s:
    :param verbosity:  (Default value = False)
    :param sanitisation_failed:  (Default value = 'error')
    :param radical:  (Default value = 'error')

    """

    assert isinstance(s, str), f"non-string smiles={s}"

    orig_smiles = s

    # if multiple molecules take the largest
    if "." in s:
        s = sorted(s.split("."), key=lambda x: len(x))[-1]

    # flatten the smiles
    stereo_smiles = s
    smiles = s.replace("@", "")
    smiles = smiles.replace("/", "")
    smiles = smiles.replace("\\", "")

    # remove isotopic stuff
    if smiles_has_isotope(smiles):
        mout.warning(f"Isotope(s) in SMILES: {smiles}")
        smiles = remove_isotopes_from_smiles(smiles)

    # canonicalise
    mol = MolFromSmiles(smiles)
    if mol:
        smiles = MolToSmiles(mol, True)
    elif sanitisation_failed == "error":
        raise SanitisationError
    elif sanitisation_failed == "warning":
        mrich.warning(f"sanitisation failed for {smiles=}")

    # check radicals
    reconstruct = False
    for i, atom in enumerate(mol.GetAtoms()):
        if not atom.GetNumRadicalElectrons():
            continue

        if radical == "warning":
            mrich.warning(f"Radical atom in {smiles=}")
        elif radical == "error":
            raise SanitisationError(f"Radical atom in {smiles=}")
        elif radical == "remove":
            mrich.warning(f"Removed radical atom")
            atom.SetNumRadicalElectrons(0)
            smiles = MolToSmiles(mol, True)
            reconstruct = True
            # atom.SetFormalCharge(0)
        else:
            raise NotImplementedError(f"Unknown option {radical=}")

    if reconstruct:
        mol = AddHs(mol)
        mol = RemoveHs(mol, implicitOnly=True)
        smiles = MolToSmiles(mol, True)
        mrich.warning(f"New {smiles=}")

    if verbosity:

        if smiles != orig_smiles:

            annotated_smiles_str = orig_smiles.replace(
                ".", f"{mcol.error}{mcol.underline}.{mcol.clear}{mcol.warning}"
            )
            annotated_smiles_str = annotated_smiles_str.replace(
                "@", f"{mcol.error}{mcol.underline}@{mcol.clear}{mcol.warning}"
            )

            mout.warning(f"SMILES was changed: {annotated_smiles_str} --> {smiles}")

    return smiles


def sanitise_mol(m):
    """

    :param m:

    """
    from rdkit.Chem import MolToMolBlock, MolFromMolBlock

    return MolFromMolBlock(MolToMolBlock(m))


def pose_gap(a, b):
    """

    :param a:
    :param b:

    """

    from numpy.linalg import norm
    from molparse.rdkit import mol_to_AtomGroup

    min_dist = None

    a = mol_to_AtomGroup(a.mol)
    b = mol_to_AtomGroup(b.mol)

    for atom1 in a.atoms:
        for atom2 in b.atoms:
            dist = norm(atom1.np_pos - atom2.np_pos)
            if min_dist is None or dist < min_dist:
                min_dist = dist

    return min_dist


ALPHANUMERIC_CHARS = "0123456789" + ascii_uppercase


def numberToBase(n, b):
    if n == 0:
        return [0]
    digits = []
    while n:
        digits.append(int(n % b))
        n //= b
    return digits[::-1]


def dt_hash():
    dt = datetime.now()
    x = int(
        dt.month * 36000 * 24 * 365.25
        + dt.day * 36000 * 24
        + dt.hour * 36000
        + dt.minute * 600
        + dt.second * 10
        + dt.microsecond / 10000
    )
    timehash = "".join([ALPHANUMERIC_CHARS[v] for v in numberToBase(x, 36)])
    return timehash


# https://github.com/rdkit/rdkit-tutorials/blob/master/notebooks/005_Chemical_space_analysis_and_visualization.ipynb

import numpy as np
from rdkit import DataStructs
from rdkit.Chem import AllChem, rdFingerprintGenerator


class FP:
    """
    Molecular fingerprint class, useful to pack features in pandas df

    Parameters
    ----------
    fp : np.array
        Features stored in numpy array
    names : list, np.array
        Names of the features
    """

    def __init__(self, fp, names):
        self.fp = fp
        self.names = names

    def __str__(self):
        return "%d bit FP" % len(self.fp)

    def __len__(self):
        return len(self.fp)


def get_cfps(
    mol, radius=2, nBits=1024, useFeatures=False, counts=False, dtype=np.float32
):
    """Calculates circular (Morgan) fingerprint.
    http://rdkit.org/docs/GettingStartedInPython.html#morgan-fingerprints-circular-fingerprints

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
    radius : float
        Fingerprint radius, default 2
    nBits : int
        Length of hashed fingerprint (without descriptors), default 1024
    useFeatures : bool
        To get feature fingerprints (FCFP) instead of normal ones (ECFP), defaults to False
    counts : bool
        If set to true it returns for each bit number of appearances of each substructure (counts). Defaults to false (fingerprint is binary)
    dtype : np.dtype
        Numpy data type for the array. Defaults to np.float32 because it is the default dtype for scikit-learn

    Returns
    -------
    ML.FP
        Fingerprint (feature) object
    """
    arr = np.zeros((1,), dtype)

    if counts is True:
        info = {}
        fp = AllChem.GetHashedMorganFingerprint(
            mol, radius, nBits, useFeatures=useFeatures
        )
        DataStructs.ConvertToNumpyArray(fp, arr)
    else:

        # https://greglandrum.github.io/rdkit-blog/posts/2023-01-18-fingerprint-generator-tutorial.html#additional-information-explaining-bits
        fmgen = rdFingerprintGenerator.GetMorganGenerator(
            radius=radius,
            fpSize=nBits,
            atomInvariantsGenerator=rdFingerprintGenerator.GetMorganFeatureAtomInvGen(),
        )

        assert not useFeatures

        DataStructs.ConvertToNumpyArray(
            fmgen.GetFingerprint(mol),
            # AllChem.GetMorganFingerprintAsBitVect(
            #     mol, radius, nBits=nBits, useFeatures=useFeatures
            # ),
            arr,
        )
    return FP(arr, range(nBits))


class SanitisationError(Exception):
    """ """

    ...
