import mrich


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
) -> "Chem.Mol":
    """Generate a header molecule for Fragalysis RHS upload"""

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
    target_access_string: str,
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
        mrich.writing(destination)
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
        "target_access_string": target_access_string,
        "trans_matrix_info": False,
    }

    mrich.print("Requesting download...")

    url = root + "/api/download_structures/"

    mrich.var("url", url)

    response = requests.post(url, json=payload)

    if response.status_code == 200:
        mrich.print("Download is ready.")
    else:
        mrich.error(f"Download request failed: {response.status_code=}")
        mrich.error(response.text)
        return None

    file_url = urllib.request.pathname2url(response.json()["file_url"])

    zip_path = destination / Path(file_url).name

    if zip_path.exists() and not overwrite:
        mrich.warning(f"Using existing {zip_path}")

    else:

        if zip_path.exists():
            mrich.warning(f"Overwriting {zip_path}")

        mrich.writing(zip_path)

        url = f"{root}/api/download_structures/?file_url={file_url}"

        mrich.var("url", url)

        filename, headers = urllib.request.urlretrieve(url, filename=zip_path)

    if not zip_path.exists():
        mrich.error("Download failed")
        return None

    if unzip:
        try:
            import zipfile

            mrich.print(f"Unzipping {zip_path}")

            target_dir = destination / Path(file_url).name.removesuffix(".zip")

            target_dir.mkdir(exist_ok=overwrite)

            mrich.writing(target_dir)
            with zipfile.ZipFile(zip_path, "r") as zip_ref:
                zip_ref.extractall(target_dir)

        except FileExistsError:
            mrich.warning(
                f"Did not unzip as directory {target_dir} exists. Set overwrite=True to override in future"
            )
            unzip = False

    mrich.success(f"Downloaded {mcol.varName}{target_name}{mcol.success} from {stack}")

    if unzip:
        return Path(target_dir)
    else:
        return Path(zip_path)


def parse_observation_longcode(longcode: str) -> dict[str]:
    """Parse a Fragalysis longcode and try to extract the following information:

    - Target name (target)
    - Crystal/dataset code (crystal)
    - Chain letter (chain)
    - Residue number (residue_number)
    - Version number (version)

    :returns: dictionary of the above keys in parentheses
    """

    import re

    dashx_split = longcode.split("-x")

    match len(dashx_split):
        case 3:

            target = dashx_split[0]

            data = dashx_split[1].removesuffix(target).split("_")

            data = [d for d in data if d]

        case 2:

            if not re.match(r"[^_]*_[A-Z]_[0-9]{3,4}_[0-9]{1,2}_.*", dashx_split[0]):
                if len(dashx_split[0]) > len(dashx_split[1]):
                    raise UnsupportedFragalysisLongcodeError(dashx_split)

                else:
                    data = dashx_split[1].split("_")[:4]
                    target = dashx_split[0]

            else:
                data = dashx_split[0].split("_")[:4]
                target = None

        case 1:

            if not re.match(
                r"[^_]*_[A-Z]_[0-9]{3,4}_[0-9]{1,2}_[^_+]*\+[A-Z]\+[0-9]{3,4}\+[0-9]{1,2}",
                dashx_split[0],
            ):
                raise UnsupportedFragalysisLongcodeError(dashx_split)

            data = dashx_split[0].split("_")[:4]
            target = None

        case _:
            mrich.var("dashx_split", dashx_split)
            raise UnsupportedFragalysisLongcodeError(dashx_split)

    match len(data):
        case 4:
            crystal, chain, residue_number, version = data
        case 3:
            crystal, chain, residue_number = data
            version = None
        case _:
            raise UnsupportedFragalysisLongcodeError(data)

    residue_number = int(residue_number)

    if len(chain) != 1:
        raise ValueError(f"{chain=}")

    version = int(version) if version else None

    return dict(
        target=target,
        crystal=crystal,
        chain=chain,
        residue_number=residue_number,
        version=version,
    )


def find_observation_longcode_matches(
    query: str, codes: list[str], debug: bool = False, allow_version_none: bool = False
) -> list[str]:

    dq = parse_observation_longcode(query)

    keys = dq.keys()

    if debug:
        mrich.var("allow_version_none", allow_version_none)
        mrich.var("dq", str(dq))

    matches = []

    for code in codes:

        if code == query:
            if debug:
                mrich.debug("exact match")
            matches.append(code)
            continue

        dc = parse_observation_longcode(code)

        for key in keys:

            if (
                allow_version_none
                and key == "version"
                and (dc[key] is None or dq[key] is None)
            ):
                continue

            if dc[key] != dq[key]:
                break
        else:
            if debug:
                mrich.debug(f"{query} matches {code}")
            matches.append(code)

    if debug:
        mrich.var("#matches", len(matches))

    if len(matches) < 1 and not allow_version_none:
        return find_observation_longcode_matches(query, codes, allow_version_none=True)

    return matches


STACK_URLS = {
    "production": "https://fragalysis.diamond.ac.uk",
    "staging": "https://fragalysis.xchem.diamond.ac.uk",
}


class UnsupportedFragalysisLongcodeError(NotImplementedError): ...
