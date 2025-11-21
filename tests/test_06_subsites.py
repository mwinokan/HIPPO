from config import *


def test_set_subsites():

    if not SUBSITES:
        return

    import hippo

    animal = hippo.HIPPO("test", DB)

    hits = animal.poses(tag="hits")
    hits.set_subsites_from_metadata_field()

    animal.db.close()


if __name__ == "__main__":
    test_set_subsites()
