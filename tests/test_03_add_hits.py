from config import *


def test_add_hits():

    if not ADD_HITS:
        return

    from pathlib import Path

    import hippo

    animal = hippo.HIPPO("test", DB)

    animal.add_hits(
        target_name=TARGET,
        aligned_directory=Path(TARGET) / "aligned_files",
        metadata_csv=Path(TARGET) / "metadata.csv",
    )

    animal.db.close()


if __name__ == "__main__":
    test_add_hits()
