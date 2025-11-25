from config import *


def test_cleanup():

    if not CLEANUP:
        return

    import os
    import shutil

    for file in CLEANUP_FILES:
        try:
            os.remove(file)
        except FileNotFoundError:
            pass

    for path in CLEANUP_DIRS:
        try:
            shutil.rmtree(path)
        except FileNotFoundError:
            pass


if __name__ == "__main__":
    test_cleanup()
