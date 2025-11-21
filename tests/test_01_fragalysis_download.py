from config import *


def test_fragalysis_download():

    if not DOWNLOAD:
        return

    from pathlib import Path
    from fragalysis.requests import download_target

    path = download_target(name=TARGET, tas=PROPOSAL)

    assert path.exists()
    assert (path / "metadata.csv").exists()
    assert (path / "aligned_files").exists()


if __name__ == "__main__":
    test_fragalysis_download()
