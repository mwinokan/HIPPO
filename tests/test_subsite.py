from config import *
from common import animal

NOT_NULL_PROPERTIES = [
    "db",
    "id",
    "table",
    "target",
    "target_id",
    "name",
    "metadata",
    "poses",
]

PROPERTIES = []


def test_properties():

    import hippo

    animal = hippo.HIPPO("test", DB)
    subsite = animal.S1

    for prop in NOT_NULL_PROPERTIES:
        value = getattr(subsite, prop)
        print(prop, value)
        assert value is not None, f"{prop} is None"

    for prop in PROPERTIES:
        value = getattr(subsite, prop)
        print(prop, value)

    animal.db.close()


if __name__ == "__main__":
    test_properties()
