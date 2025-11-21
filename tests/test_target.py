from config import *
from common import animal

NOT_NULL_PROPERTIES = [
    "id",
    "name",
    "feature_ids",
    "features",
    "subsites",
]

PROPERTIES = []


def test_properties():

    import hippo

    animal = hippo.HIPPO("test", DB)
    target = animal.T1

    for prop in NOT_NULL_PROPERTIES:
        value = getattr(target, prop)
        print(prop, value)
        assert value is not None, f"{prop} is None"

    for prop in PROPERTIES:
        value = getattr(target, prop)
        print(prop, value)

    animal.db.close()


if __name__ == "__main__":
    test_properties()
