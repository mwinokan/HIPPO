import logging


logger = logging.getLogger("HIPPO")


def generate_header(
    pose,
    method,
    ref_url,
    submitter_name,
    submitter_email,
    submitter_institution,
    generation_date: str | None = None,
    extras=None,
    metadata: bool = True,
):

    extras = extras or {}

    from rdkit.Chem.AllChem import EmbedMolecule
    from molparse.rdkit import mol_from_smiles
    from datetime import date

    header = mol_from_smiles(pose.compound.smiles)

    header.SetProp("_Name", "ver_1.2")
    EmbedMolecule(header)

    generation_date = str(generation_date or date.today())

    header.SetProp("ref_url", ref_url)
    header.SetProp("submitter_name", submitter_name)
    header.SetProp("submitter_email", submitter_email)
    header.SetProp("submitter_institution", submitter_institution)
    header.SetProp("generation_date", generation_date)
    header.SetProp("method", method)

    if metadata:
        for k, v in pose.metadata.items():
            header.SetProp(k, str(k))

    for k, v in extras.items():
        header.SetProp(k, str(v))

    return header


def download_target(
    target_name: str,
    *,
    destination: "str | Path" = ".",
    stack: str = "production",
    unzip: bool = True,
    overwrite: bool = False,
) -> "Path":
    """Download a target from Fragalysis

    :param target_name: Name of the target
    :param destination: where to put the file(s)
    :param stack: Choose the Fragalysis stack from ['production', 'staging']. Defaults to 'production'
    :param unzip: Unzip the download
    :param overwrite: Overwrite existing downloads
    :returns: a pathlib.Path object to the .zip archive or unpacked directory

    """

    import mcol
    import requests
    from pathlib import Path
    import urllib.request

    destination = Path(destination)

    if not destination.exists():
        logger.writing(destination)
        destination.mkdir(exist_ok=True, parents=True)

    root = STACK_URLS[stack]

    payload = {
        "all_aligned_structures": True,
        "cif_info": False,
        "diff_file": False,
        "event_file": False,
        "file_url": "",
        "map_info": False,
        "metadata_info": True,
        "mtz_info": False,
        "pdb_info": False,
        "proteins": "",
        "sigmaa_file": False,
        "single_sdf_file": True,
        "static_link": False,
        "target_name": target_name,
        "trans_matrix_info": False,
    }

    logger.info("Requesting download...")

    url = root + "/api/download_structures/"

    logger.var("url", url)

    response = requests.post(url, json=payload)

    if response.status_code == 200:
        logger.info("Download is ready.")
    else:
        logger.error(f"Download request failed: {response.status_code=}")
        logger.error(response.text)
        return None

    file_url = urllib.request.pathname2url(response.json()["file_url"])

    zip_path = destination / Path(file_url).name

    if zip_path.exists() and not overwrite:
        logger.warning(f"Using existing {zip_path}")

    else:

        if zip_path.exists():
            logger.warning(f"Overwriting {zip_path}")

        logger.writing(zip_path)

        url = f"{root}/api/download_structures/?file_url={file_url}"

        logger.var("url", url)

        filename, headers = urllib.request.urlretrieve(url, filename=zip_path)

    if not zip_path.exists():
        logger.error("Download failed")
        return None

    if unzip:
        try:
            import zipfile

            logger.info(f"Unzipping {zip_path}")

            target_dir = destination / Path(file_url).name.removesuffix(".zip")

            target_dir.mkdir(exist_ok=overwrite)

            logger.writing(target_dir)
            with zipfile.ZipFile(zip_path, "r") as zip_ref:
                zip_ref.extractall(target_dir)

        except FileExistsError:
            logger.warning(
                f"Did not unzip as directory {target_dir} exists. Set overwrite=True to override in future"
            )
            unzip = False

    logger.success(f"Downloaded {mcol.varName}{target_name}{mcol.success} from {stack}")

    if unzip:
        return Path(target_dir)
    else:
        return Path(zip_path)


STACK_URLS = {
    "production": "https://fragalysis.diamond.ac.uk",
    "staging": "https://fragalysis.xchem.diamond.ac.uk",
}
