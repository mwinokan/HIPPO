from config import *


def test_setup_animal():

    if not SETUP:
        return

    from pathlib import Path
    import hippo

    animal = hippo.HIPPO("test", DB)
    animal.summary()

    if isinstance(DB, str):
        assert Path(DB).exists()

    animal.db.close()


if __name__ == "__main__":
    test_setup_animal()
