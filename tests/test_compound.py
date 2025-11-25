from config import *

NOT_NULL_PROPERTIES = [
    "id",
    "inchikey",
    "name",
    "smiles",
    "mol",
    "num_heavy_atoms",
    "molecular_weight",
    "num_rings",
    "formula",
    "atomtype_dict",
    "metadata",
    "db",
    "tags",
    "poses",
    "best_placed_pose",
    "num_poses",
    "num_reactions",
    "num_reactant",
    "num_scaffolds",
    "dict",
    "is_scaffold",
    "is_elab",
    "is_product",
    "table",
]

PROPERTIES = [
    "alias",
    "elabs",
    "reaction",
    "reactions",
    "scaffolds",
    "num_atoms_added",
]


def test_properties():

    import hippo

    animal = hippo.HIPPO("test", DB)
    compound = animal.C1

    for prop in NOT_NULL_PROPERTIES:
        value = getattr(compound, prop)
        print(prop, value)
        assert value is not None, f"{prop} is None"

    for prop in PROPERTIES:
        value = getattr(compound, prop)
        print(prop, value)

    animal.db.close()


if __name__ == "__main__":
    test_properties()
