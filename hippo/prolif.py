"""Functions for ProLIF interaction profiling"""

import mrich
from mrich import print

import molparse as mp

from molparse.rdkit.features import FEATURE_FAMILIES, INTERACTION_TYPES

INTERACTION_TYPES = list(INTERACTION_TYPES.values()) + ["VdWContact"]
FEATURE_FAMILIES = list(FEATURE_FAMILIES) + ["VdWSphere"]


def guess_feature_families(
    interaction_type: str, prot_atom_names: list[str], lig_atom_ids: list[int]
) -> tuple[str, str, str]:
    """Guess feature families from atom names

    :param interaction_type: interaction type

    """

    match interaction_type:
        case "VdWContact":
            lig_feature_family = "VdWSphere"
            prot_feature_family = "VdWSphere"

        case "Hydrophobic":
            lig_feature_family = (
                "LumpedHydrophobe" if len(lig_atom_ids) > 1 else "Hydrophobe"
            )
            prot_feature_family = (
                "LumpedHydrophobe" if len(prot_atom_names) > 1 else "Hydrophobe"
            )

        case "Anionic":
            lig_feature_family = "NegIonizable"
            prot_feature_family = "PosIonizable"
            interaction_type = "Electrostatic"

        case "Cationic":
            lig_feature_family = "PosIonizable"
            prot_feature_family = "NegIonizable"
            interaction_type = "Electrostatic"

        case "CationPi":
            lig_feature_family = "PosIonizable"
            prot_feature_family = "Aromatic"
            interaction_type = "π-cation"

        case "PiCation":
            lig_feature_family = "Aromatic"
            prot_feature_family = "PosIonizable"
            interaction_type = "π-cation"

        case "PiStacking":
            lig_feature_family = "Aromatic"
            prot_feature_family = "Aromatic"
            interaction_type = "π-stacking"

        case "PiStacking":
            lig_feature_family = "Aromatic"
            prot_feature_family = "Aromatic"
            interaction_type = "π-stacking"

        case "EdgeToFace":
            lig_feature_family = "Aromatic"
            prot_feature_family = "Aromatic"
            interaction_type = "π-stacking (EdgeToFace)"

        case "FaceToFace":
            lig_feature_family = "Aromatic"
            prot_feature_family = "Aromatic"
            interaction_type = "π-stacking (FaceToFace)"

        case "HBAcceptor":
            lig_feature_family = "Acceptor"
            prot_feature_family = "Donor"
            interaction_type = "Hydrogen Bond"

        case "HBDonor":
            lig_feature_family = "Donor"
            prot_feature_family = "Acceptor"
            interaction_type = "Hydrogen Bond"

        case "XBAcceptor":
            lig_feature_family = "Acceptor"
            prot_feature_family = "Donor"
            interaction_type = "Halogen Bond"

        case "XBDonor":
            lig_feature_family = "Donor"
            prot_feature_family = "Acceptor"
            interaction_type = "Halogen Bond"

        case "MetalAcceptor":
            lig_feature_family = None
            prot_feature_family = "Metal"
            interaction_type = "Metal complexation"

        case "MetalDonor":
            lig_feature_family = "Metal"
            prot_feature_family = None
            interaction_type = "Metal complexation"

        case _:
            prot_feature_family = None
            lig_feature_family = None

            mrich.error("Unsupported interaction_type", interaction_type)

    assert interaction_type in INTERACTION_TYPES
    assert prot_feature_family in FEATURE_FAMILIES or prot_feature_family is None
    assert lig_feature_family in FEATURE_FAMILIES or lig_feature_family is None

    return prot_feature_family, lig_feature_family, interaction_type


def parse_prolif_interactions(
    pose: "Pose",
    fp: "plf.Fingerprint",
    protonated_sys: "molparse.System",
    table: str = "temp_interaction",
    debug: bool = False,
) -> None:
    """Parse ProLIF output into HIPPO database table"""

    target_id = pose.target.id

    for key, value in fp.ifp[0].items():

        res_number = key[1].number
        res_name = key[1].name
        chain_name = key[1].chain

        chain = protonated_sys.get_chain(chain_name)
        residue = chain.residues.get_matches(number=res_number)

        prot_mol = residue.rdkit_mol
        prot_group = mp.AtomGroup.from_pdb_block(mp.rdkit.mol_to_pdb_block(prot_mol))

        if not residue.name == res_name:
            mrich.debug(protonated_sys.name + ".pdb")
            raise AssertionError(
                f"Residue name mismatch: [sys]={residue.name} [prolif]={res_name}"
            )

        for interaction_type, interaction_dicts in value.items():

            for interaction_dict in interaction_dicts:

                angle = interaction_dict.get("angle")
                distance = interaction_dict.get("distance")

                lig_atom_ids = list(interaction_dict["indices"]["ligand"])

                # insert a dummy protein feature
                prot_interaction_atoms = [
                    prot_group.atoms[i] for i in interaction_dict["indices"]["protein"]
                ]
                prot_atom_names = [a.name for a in prot_interaction_atoms]

                # rename stuff
                prot_family, lig_family, interaction_type = guess_feature_families(
                    interaction_type, prot_atom_names, lig_atom_ids
                )

                feature_id = pose.db.insert_feature(
                    family=prot_family,
                    target=target_id,
                    chain_name=chain_name,
                    residue_name=res_name,
                    residue_number=residue.number,
                    atom_names=prot_atom_names,
                    commit=False,
                )

                if not feature_id:

                    sql = f"""
                    feature_target = {target_id} 
                    AND feature_family = '{prot_family}' 
                    AND feature_chain_name = '{chain_name}' 
                    AND feature_residue_name = '{res_name}' 
                    AND feature_residue_number = {residue.number}
                    AND feature_atom_names = '{" ".join(sorted(prot_atom_names))}'
                    """

                    try:
                        (feature_id,) = pose.db.select_id_where(
                            table="feature", key=sql
                        )
                    except:
                        feature_id = pose.db.insert_feature(
                            family=prot_family,
                            target=target_id,
                            chain_name=chain_name,
                            residue_name=res_name,
                            residue_number=residue.number,
                            atom_names=prot_atom_names,
                            warn_duplicate=True,
                            commit=False,
                        )
                        raise

                # insert into the Database
                interaction_id = pose.db.insert_interaction(
                    feature=feature_id,
                    pose=pose.id,
                    type=interaction_type,
                    family=lig_family,
                    atom_ids=[i + 1 for i in lig_atom_ids],
                    prot_coord=None,
                    lig_coord=None,
                    distance=distance,
                    angle=angle,
                    energy=None,
                    commit=False,
                    table=table,
                )

                if debug:
                    mrich.debug(
                        f"Residue: {res_name} {res_number} {chain_name}. Interaction: {interaction_type}. Ligand: {lig_family} ({lig_atom_ids}). Protein: {prot_family} ({prot_atom_names})"
                    )
