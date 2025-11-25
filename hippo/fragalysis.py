"""Functions for interfacing with Fragalysis data"""

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

    match = re.search(
        r"(.*)_([A-z]_[0-9]*_[0-9])_(.*)\+([A-z]\+[0-9]*\+[0-9])_.LIG", longcode
    )

    if not match:
        raise UnsupportedFragalysisLongcodeError(longcode)

    cryst_str, lig_str, _, _ = match.groups()

    chain, residue_number, version = lig_str.split("_")

    residue_number = int(residue_number)
    version = int(version)

    if match := re.search(r"(.*)-(\w[0-9]{4})", cryst_str):

        target_name = match.group(0)
        crystal = match.group(1)

    else:

        target_name = None
        crystal = cryst_str

    return dict(
        target=target_name,
        crystal=crystal,
        chain=chain,
        residue_number=residue_number,
        version=version,
    )


def find_observation_longcode_matches(
    query: str, codes: list[str], debug: bool = False, allow_version_none: bool = False
) -> list[str]:
    """find_observation_longcode_matches"""

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


class UnsupportedFragalysisLongcodeError(NotImplementedError):
    """Provided Fragalysis observation long code syntax is not supported"""

    ...
