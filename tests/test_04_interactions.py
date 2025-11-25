from config import *


def test_calculate_interactions():

    import hippo

    animal = hippo.HIPPO("test", DB)

    for pose in animal.poses:
        pose.calculate_interactions()

    animal.db.close()


if __name__ == "__main__":
    test_calculate_interactions()
