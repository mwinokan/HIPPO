from config import *
import hippo


def test_calculate_all_scaffolds():
    if SCAFFOLDS:
        animal = hippo.HIPPO("test", DB)
        animal.db.calculate_all_scaffolds()
        animal.db.close()


def test_calculate_all_murcko_scaffolds():
    if SCAFFOLDS:
        animal = hippo.HIPPO("test", DB)
        animal.db.calculate_all_murcko_scaffolds()
        animal.db.close()


if __name__ == "__main__":
    test_calculate_all_scaffolds()
    test_calculate_all_murcko_scaffolds()
