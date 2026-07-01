from config import *

NOT_NULL_PROPERTIES = [
    'unique',
]

PROPERTIES = []


def test_properties():

    import hippo

    animal = hippo.HIPPO('test', DB)
    tag_table = animal.tags

    for prop in NOT_NULL_PROPERTIES:
        value = getattr(tag_table, prop)
        print(prop, value)
        assert value is not None, f'{prop} is None'

    for prop in PROPERTIES:
        value = getattr(tag_table, prop)
        print(prop, value)

    animal.db.close()


def test_summary():

    import hippo

    animal = hippo.HIPPO('test', DB)
    tag_table = animal.tags
    tag_table.summary()


if __name__ == '__main__':
    test_properties()
    test_summary()
