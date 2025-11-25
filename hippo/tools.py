"""Generic tools for use in the HIPPO package"""

import re
import numpy as np
from molparse.rdkit import mol_from_smiles
from rdkit.Chem.inchi import MolToInchiKey
from rdkit.Chem import MolFromSmiles, MolToSmiles, AddHs, RemoveHs
import mcol
from datetime import datetime
from string import ascii_uppercase

import mrich


def df_row_to_dict(df_row) -> dict:
    """Convert a dataframe row to a dictionary

    :param df_row: pandas dataframe row / series
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


def remove_other_ligands(
    sys: "molparse.System", residue_number: int, chain: str
) -> "molparse.System":
    """Remove ligands other than the specified one"""

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


def inchikey_from_smiles(smiles: str) -> str:
    """InChI-Key from smiles string"""
    mol = mol_from_smiles(smiles)
    return MolToInchiKey(mol)


def flat_inchikey(smiles: str) -> str:
    """Stereochemistry-flattened InChI-Key from smiles string"""
    smiles = sanitise_smiles(smiles)
    return inchikey_from_smiles(smiles)


def remove_isotopes_from_smiles(smiles: str) -> str:
    """Remove isotopes from smiles string"""

    mol = MolFromSmiles(smiles)

    atom_data = [(atom, atom.GetIsotope()) for atom in mol.GetAtoms()]

    for atom, isotope in atom_data:
        if isotope:
            atom.SetIsotope(0)

    return MolToSmiles(mol)


def smiles_has_isotope(smiles: str, regex=True) -> bool:
    """Does provided smiles string contain isotopes?"""
    if regex:
        return re.search(r"([\[][0-9]+[A-Z]+\])", smiles)
    else:
        mol = MolFromSmiles(smiles)
        return any(atom.GetIsotope() for atom in mol.GetAtoms())


REPLACE = {
    "[STB]": "[S]",
}


def sanitise_smiles(
    s: str,
    verbosity: bool = False,
    sanitisation_failed: str = "error",
    radical: str = "error",
) -> str:
    """Sanitise smiles by:

    - Taking largest fragment
    - Flattening stereochemistry
    - Removing isotopes
    - RDKit round-trip
    - Treating radicals

    :param s: input smiles string
    :param verbosity: print smiles changes (Default value = False)
    :param sanitisation_failed: behvaiour when sanitisation fails, choose from ["error", "warning", "quiet"] (Default value = 'error')
    :param radical: behvaiour when radicals occur, choose from ["error", "warning", "remove"] (Default value = 'error')
    :returns: SMILES string
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
        mrich.warning(f"Isotope(s) in SMILES: {smiles}")
        smiles = remove_isotopes_from_smiles(smiles)

    # replace specific sequences
    for key in REPLACE:
        if key in smiles:
            smiles = smiles.replace(key, REPLACE[key])

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
    for atom in mol.GetAtoms():
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

            mrich.warning(f"SMILES was changed: {annotated_smiles_str} --> {smiles}")

    return smiles


def sanitise_mol(m: "rdkit.Chem.Mol") -> "rdkit.Chem.Mol":
    """Sanitise by RDKit round-trip"""
    from rdkit.Chem import MolToMolBlock, MolFromMolBlock

    return MolFromMolBlock(MolToMolBlock(m))


def pose_gap(a: "Pose", b: "Pose") -> float:
    """Calculate minimum distance between two :class:`.Pose` objects"""

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


def number_to_base(n: int, b: int) -> int:
    """Convert an integer `n` into base `b` representation"""
    if n == 0:
        return [0]
    digits = []
    while n:
        digits.append(int(n % b))
        n //= b
    return digits[::-1]


def dt_hash() -> str:
    """Create 7 alphanumeric-character hash based on current timestamp"""
    dt = datetime.now()
    x = int(
        dt.month * 36000 * 24 * 365.25
        + dt.day * 36000 * 24
        + dt.hour * 36000
        + dt.minute * 600
        + dt.second * 10
        + dt.microsecond / 10000
    )
    timehash = "".join([ALPHANUMERIC_CHARS[v] for v in number_to_base(x, 36)])
    return f"{timehash:>07}"


class SanitisationError(Exception):
    """Something went wrong in Molecule/SMILES sanitisation"""

    ...
