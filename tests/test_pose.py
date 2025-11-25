from config import *
from common import animal

NOT_NULL_PROPERTIES = [
    "db",
    "id",
    "inchikey",
    "alias",
    "name",
    "smiles",
    "target",
    "compound_id",
    "compound",
    "path",
    "mol",
    "protonated_mol",
    "protein_system",
    "complex_system",
    "has_complex_pdb_path",
    "metadata",
    "has_fingerprint",
    "tags",
    "features",
    "dict",
    "table",
    "num_heavy_atoms",
    "num_scaffolds",
    "scaffold_ids",
    "interactions",
    "classic_fingerprint",
    "mol_path",
    "apo_path",
]

PROPERTIES = [
    "reference",
    "reference_id",
    "inspirations",
    "derivatives",
    "num_atoms_added",
    "num_atoms_added_wrt_scaffolds",
    "num_atoms_added_wrt_inspirations",
    "energy_score",
    "distance_score",
    "inspiration_score",
    "subsites",
]


def test_properties():

    import hippo

    animal = hippo.HIPPO("test", DB)
    pose = animal.P1

    for prop in NOT_NULL_PROPERTIES:
        value = getattr(pose, prop)
        print(prop, value)
        assert value is not None, f"{prop} is None"

    for prop in PROPERTIES:
        value = getattr(pose, prop)
        print(prop, value)
    animal.db.close()


if __name__ == "__main__":
    test_properties()
