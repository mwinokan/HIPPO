"""Functions for interfacing with Fragalysis data"""

import os
import re
from dataclasses import dataclass, fields

import mrich
from rdkit import Chem

GENERATED_TAG_COLS = [
    'ConformerSites alias',
    'CanonSites alias',
    'CrystalformSites alias',
    'Quatassemblies alias',
    'Crystalforms alias',
    'ConformerSites upload name',
    'CanonSites upload name',
    'CrystalformSites upload name',
    'Quatassemblies upload name',
    'Crystalforms upload name',
    'ConformerSites short tag',
    'CanonSites short tag',
    'CrystalformSites short tag',
    'Quatassemblies short tag',
    'Crystalforms short tag',
    'Centroid res',
    'Experiment code',
    'PoseModel',
]


META_IGNORE_COLS = [
    'Code',
    'Long code',
    'CompoundModel code',
    'Smiles',
    'Downloaded',
    'Main status',
    'GOOD count',
    'MEDIOCRE count',
    'BAD count',
    'RefinementResolution',
]


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
) -> Chem.rdchem.Mol:
    """Generate a header molecule for Fragalysis RHS upload"""

    extras = extras or {}

    from datetime import date

    from molparse.rdkit import mol_from_smiles
    from rdkit.Chem.AllChem import EmbedMolecule

    header = mol_from_smiles(pose.compound.compound_smiles)

    header.SetProp('_Name', 'ver_1.2')
    EmbedMolecule(header)

    generation_date = str(generation_date or date.today())

    header.SetProp('ref_url', ref_url)
    header.SetProp('submitter_name', submitter_name)
    header.SetProp('submitter_email', submitter_email)
    header.SetProp('submitter_institution', submitter_institution)
    header.SetProp('generation_date', generation_date)
    header.SetProp('method', method)

    if metadata:
        for k, _ in pose.pose_metadata.items():
            header.SetProp(k, str(k))

    for k, v in extras.items():
        header.SetProp(k, str(v))

    return header


@dataclass
class LongcodeRecord:
    crystal: str  # full crystal token, e.g. A71EV2A-x0152
    protein_name: str  # e.g. A71EV2A
    chain: str
    residue_number: int
    version: int
    altloc: str | None = None


# Current Fragalysis longcode format. The code combines two sites
# (e.g. A71EV2A-x0152_A_147_0_1_A71EV2A-x0526+A+147+0+1__LIG); we parse only the
# first (underscore-separated) site and ignore the second (plus-separated) one.
# First-site groups: crystal token, then chain_resnum_altloc_version.
_LONGCODE_RE = re.compile(
    r'(.*)_'  # crystal token, e.g. A71EV2A-x0152
    r'([A-Za-z]+_[0-9]+_[A-Za-z0-9]+_[0-9]+)_'  # chain_resnum_altloc_version
    r'.*\+[A-Za-z]+\+[0-9]+\+[A-Za-z0-9]+\+[0-9]+'  # second site (ignored)
    r'_.LIG'
)

# DEPRECATED(longcode-altloc): pre-altloc Fragalysis longcode format
# (e.g. D68EV3CPROB-x0455_A_209_1_7gp9+A+201+1__LIG, no altloc group). Remove this
# regex and its branch in parse_observation_longcode once all data uses the
# current format above.
_LONGCODE_RE_LEGACY = re.compile(
    r'(.*)_([A-Za-z]_[0-9]*_[0-9])_(.*)\+([A-Za-z]\+[0-9]*\+[0-9])_.LIG'
)

# Split a crystal token (e.g. A71EV2A-x0152) into protein name and crystal id.
_CRYSTAL_RE = re.compile(r'(.*)-(\w[0-9]{4})')


def parse_observation_longcode(longcode: str) -> LongcodeRecord:
    """Parse a Fragalysis observation longcode (first site only).

    Extracts:

    - Crystal token (crystal), e.g. ``A71EV2A-x0152``
    - Protein name (protein_name), e.g. ``A71EV2A``
    - Chain letter (chain)
    - Residue number (residue_number)
    - Altloc (altloc; ``None`` for the older format)
    - Version number (version)
    """

    altloc = None

    match = _LONGCODE_RE.search(longcode)
    if match:
        cryst_str = match.group(1)
        chain, residue_number, altloc, version = match.group(2).split('_')
    else:
        # DEPRECATED(longcode-altloc): pre-altloc format had no altloc group;
        # remove this branch (and _LONGCODE_RE_LEGACY) once all data uses the
        # current format.
        match = _LONGCODE_RE_LEGACY.search(longcode)
        if not match:
            raise UnsupportedFragalysisLongcodeError(longcode)
        cryst_str = match.group(1)
        chain, residue_number, version = match.group(2).split('_')
        # end DEPRECATED(longcode-altloc)

    residue_number = int(residue_number)
    version = int(version)

    if m := _CRYSTAL_RE.search(cryst_str):
        crystal = m.group(0)  # full token, e.g. A71EV2A-x0152
        protein_name = m.group(1)  # e.g. A71EV2A
    else:
        crystal = cryst_str
        protein_name = ''

    return LongcodeRecord(
        crystal=crystal,
        protein_name=protein_name,
        chain=chain,
        residue_number=residue_number,
        version=version,
        altloc=altloc,
    )


def find_observation_longcode_matches(
    query: str, codes: list[str], debug: bool = False, allow_version_none: bool = False
) -> list[str]:
    """find_observation_longcode_matches"""

    dq = parse_observation_longcode(query)

    if debug:
        mrich.var('allow_version_none', allow_version_none)
        mrich.var('dq', str(dq))

    matches = []

    for code in codes:
        if code == query:
            if debug:
                mrich.debug('exact match')
            matches.append(code)
            continue

        dc = parse_observation_longcode(code)

        for key in fields(dq):
            if (
                allow_version_none
                and key.name == 'version'
                and (getattr(dc, key.name) is None or getattr(dq, key.name) is None)
            ):
                continue

            if getattr(dc, key.name) != getattr(dq, key.name):
                break
        else:
            if debug:
                mrich.debug(f'{query} matches {code}')
            matches.append(code)

    if debug:
        mrich.var('#matches', len(matches))

    if len(matches) < 1 and not allow_version_none:
        return find_observation_longcode_matches(query, codes, allow_version_none=True)

    return matches


STACK_URLS = {
    'production': 'https://fragalysis.diamond.ac.uk',
    'staging': 'https://fragalysis.xchem.diamond.ac.uk',
    # testing
    'localhost': 'http://localhost:8080',
}

# Developers can add or override stacks via the environment without editing this
# shared file. Each ``HIPPO_STACK_URL_<NAME>`` variable becomes the ``<name>``
# stack (lower-cased), e.g. set in your (gitignored) .env:
#     HIPPO_STACK_URL_DOCKERHOST=http://host.docker.internal:8080
# then use stack='dockerhost'.
_STACK_URL_ENV_PREFIX = 'HIPPO_STACK_URL_'
STACK_URLS.update(
    {
        key[len(_STACK_URL_ENV_PREFIX) :].lower(): value
        for key, value in os.environ.items()
        if key.startswith(_STACK_URL_ENV_PREFIX) and value
    }
)


class UnsupportedFragalysisLongcodeError(NotImplementedError):
    """Provided Fragalysis observation long code syntax is not supported"""

    ...
