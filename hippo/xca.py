import mrich
from mrich import print

import re


def parse_observation_longcode(longcode: str) -> dict[str]:
    """Parse an XChemAlign longcode and try to extract the following information:

    - Target name (target)
    - Crystal/dataset code (crystal)
    - Chain letter (chain)
    - Residue number (residue_number)
    - Version number (version)

    :returns: dictionary of the above keys in parentheses
    """

    match = re.search(
        r"^(.*)-(.\d{4})_(.)_(\d*)_(\d)_.*-.\d{4}\+.\+\d*\+\d_.LIG$", longcode
    )

    if not match:
        raise UnsupportedXCALongcodeError(longcode)

    target_name, crystal, chain, residue_number, version = match.groups()

    return dict(
        target=target_name,
        crystal=crystal,
        chain=chain,
        residue_number=int(residue_number),
        version=int(version),
    )


class UnsupportedXCALongcodeError(NotImplementedError): ...
