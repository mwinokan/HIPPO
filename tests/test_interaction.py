from config import *
from common import animal

NOT_NULL_PROPERTIES = [
    "id",
    "table",
    "db",
    "family",
    "pose_id",
    "pose",
    "feature_id",
    "feature",
    "residue_name",
    "residue_number",
    "atom_ids",
    "prot_coord",
    "lig_coord",
    "distance",
    "family_str",
    "type",
    "description",
]

PROPERTIES = [
    "angle",
    "energy",
]


def test_properties():

    import hippo

    animal = hippo.HIPPO("test", DB)
    interaction = animal.I1

    for prop in NOT_NULL_PROPERTIES:
        value = getattr(interaction, prop)
        print(prop, value)
        assert value is not None, f"{prop} is None"

    for prop in PROPERTIES:
        value = getattr(interaction, prop)
        print(prop, value)

    animal.db.close()


if __name__ == "__main__":
    test_properties()
