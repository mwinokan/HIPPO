from config import *
from common import animal

NOT_NULL_PROPERTIES = [
    "id",
    "family",
    "target",
    "chain_name",
    "residue_name",
    "residue_number",
    "atom_names",
]

PROPERTIES = []


def test_properties():

    import hippo

    animal = hippo.HIPPO("test", DB)
    feature = animal.F1

    for prop in NOT_NULL_PROPERTIES:
        value = getattr(feature, prop)
        print(prop, value)
        assert value is not None, f"{prop} is None"

    for prop in PROPERTIES:
        value = getattr(feature, prop)
        print(prop, value)

    animal.db.close()


if __name__ == "__main__":
    test_properties()
