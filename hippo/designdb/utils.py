"""Generic tools for use in the HIPPO package"""

import ast
import cProfile
import functools
import json
import os
import re
from collections.abc import Iterator
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
from datetime import datetime
from string import ascii_uppercase
from typing import TYPE_CHECKING

import mcol
import molparse as mp
import mrich
import numpy as np
from django.db import ProgrammingError, connection
from django.db.models import Aggregate, OuterRef, Subquery
from molparse.rdkit import mol_from_smiles
from rdkit import Chem
from rdkit.Chem import AddHs, MolFromSmiles, MolToSmiles, RegistrationHash, RemoveHs
from rdkit.Chem.inchi import MolToInchiKey
from rdkit.Chem.MolStandardize import rdMolStandardize

if TYPE_CHECKING:
    from designdb.models import PoseModel


def strip_sql(sql) -> str:
    """Reduce unecessary whitespace in SQL"""
    return re.sub(r'\s+', ' ', sql).strip()


def df_row_to_dict(df_row) -> dict:
    """Convert a dataframe row to a dictionary

    :param df_row: pandas dataframe row / series
    """

    assert len(df_row) == 1, f'{len(df_row)=}'

    data = {}

    for col in df_row.columns:
        if col == 'Unnamed: 0':
            continue

        value = df_row[col].values[0]

        if not isinstance(value, str) and np.isnan(value):
            value = None

        data[col] = value

    return data


def remove_other_ligands(sys: mp.System, residue_number: int, chain: str) -> mp.System:
    """Remove ligands other than the specified one"""

    ligand_residues = [r.number for r in sys['rLIG'] if r.number != residue_number]

    # if ligand_residues:
    for c in sys.chains:
        if c.name != chain:
            c.remove_residues(names=['LIG'], verbosity=0)
        elif ligand_residues:
            c.remove_residues(numbers=ligand_residues, verbosity=0)

    # print([r.name_number_str for r in sys['rLIG']])

    assert len([r.name_number_str for r in sys['rLIG']]) == 1, (
        f'{sys.name} {[r.name_number_str for r in sys["rLIG"]]}'
    )

    return sys


def inchikey_from_smiles(smiles: str) -> str:
    """InChI-Key from smiles string"""
    mol = mol_from_smiles(smiles)
    return MolToInchiKey(mol)


def flat_inchikey(smiles: str) -> str:
    """Stereochemistry-flattened InChI-Key from smiles string"""
    smiles = sanitise_smiles(smiles)
    return inchikey_from_smiles(smiles)


def remove_isotopes_from_smiles(smiles: str) -> str:
    """Remove isotopes from smiles string"""

    mol = MolFromSmiles(smiles)

    atom_data = [(atom, atom.GetIsotope()) for atom in mol.GetAtoms()]

    for atom, isotope in atom_data:
        if isotope:
            atom.SetIsotope(0)

    return MolToSmiles(mol)


def smiles_has_isotope(smiles: str, regex=True) -> bool:
    """Does provided smiles string contain isotopes?"""
    if regex:
        return re.search(r'([\[][0-9]+[A-Z]+\])', smiles)
    else:
        mol = MolFromSmiles(smiles)
        return any(atom.GetIsotope() for atom in mol.GetAtoms())


REPLACE = {
    '[STB]': '[S]',
}


def sanitise_smiles(
    s: str,
    verbosity: bool = False,
    sanitisation_failed: str = 'error',
    radical: str = 'error',
) -> str:
    """Sanitise smiles by:

    - Taking largest fragment
    - Flattening stereochemistry
    - Removing isotopes
    - RDKit round-trip
    - Treating radicals

    :param s: input smiles string
    :param verbosity: print smiles changes (Default value = False)
    :param sanitisation_failed: behvaiour when sanitisation fails,
        choose from ["error", "warning", "quiet"] (Default value = 'error')
    :param radical: behvaiour when radicals occur, choose from
        ["error", "warning", "remove"] (Default value = 'error')
    :returns: SMILES string
    """

    assert isinstance(s, str), f'non-string smiles={s}'

    orig_smiles = s

    # if multiple molecules take the largest
    if '.' in s:
        s = sorted(s.split('.'), key=lambda x: len(x))[-1]

    # flatten the smiles
    # stereo_smiles = s
    smiles = s.replace('@', '')
    smiles = smiles.replace('/', '')
    smiles = smiles.replace('\\', '')

    # remove isotopic stuff
    if smiles_has_isotope(smiles):
        mrich.warning(f'Isotope(s) in SMILES: {smiles}')
        smiles = remove_isotopes_from_smiles(smiles)

    # replace specific sequences
    for key in REPLACE:
        if key in smiles:
            smiles = smiles.replace(key, REPLACE[key])

    # canonicalise
    mol = MolFromSmiles(smiles)
    if mol:
        smiles = MolToSmiles(mol, True)
    elif sanitisation_failed == 'error':
        raise SanitisationError
    elif sanitisation_failed == 'warning':
        mrich.warning(f'sanitisation failed for {smiles=}')

    # check radicals
    reconstruct = False
    for atom in mol.GetAtoms():
        if not atom.GetNumRadicalElectrons():
            continue

        if radical == 'warning':
            mrich.warning(f'Radical atom in {smiles=}')
        elif radical == 'error':
            raise SanitisationError(f'Radical atom in {smiles=}')
        elif radical == 'remove':
            mrich.warning('Removed radical atom')
            atom.SetNumRadicalElectrons(0)
            smiles = MolToSmiles(mol, True)
            reconstruct = True
            # atom.SetFormalCharge(0)
        else:
            raise NotImplementedError(f'Unknown option {radical=}')

    if reconstruct:
        mol = AddHs(mol)
        mol = RemoveHs(mol, implicitOnly=True)
        smiles = MolToSmiles(mol, True)
        mrich.warning(f'New {smiles=}')

    if verbosity:
        if smiles != orig_smiles:
            annotated_smiles_str = orig_smiles.replace(
                '.', f'{mcol.error}{mcol.underline}.{mcol.clear}{mcol.warning}'
            )
            annotated_smiles_str = annotated_smiles_str.replace(
                '@', f'{mcol.error}{mcol.underline}@{mcol.clear}{mcol.warning}'
            )

            mrich.warning(f'SMILES was changed: {annotated_smiles_str} --> {smiles}')

    return smiles


def sanitise_mol(m: Chem.rdchem.Mol) -> Chem.rdchem.Mol:
    """Sanitise by RDKit round-trip"""
    from rdkit.Chem import MolFromMolBlock, MolToMolBlock

    return MolFromMolBlock(MolToMolBlock(m))


def pose_gap(a: 'PoseModel', b: 'PoseModel') -> float:
    """Calculate minimum distance between two :class:`.PoseModel` objects"""

    from molparse.rdkit import mol_to_AtomGroup
    from numpy.linalg import norm

    # avoiding circular imports

    min_dist = None

    a = mol_to_AtomGroup(a.mol)
    b = mol_to_AtomGroup(b.mol)

    for atom1 in a.atoms:
        for atom2 in b.atoms:
            dist = norm(atom1.np_pos - atom2.np_pos)
            if min_dist is None or dist < min_dist:
                min_dist = dist

    return min_dist


ALPHANUMERIC_CHARS = '0123456789' + ascii_uppercase


def number_to_base(n: int, b: int) -> int:
    """Convert an integer `n` into base `b` representation"""
    if n == 0:
        return [0]
    digits = []
    while n:
        digits.append(int(n % b))
        n //= b
    return digits[::-1]


def dt_hash() -> str:
    """Create 7 alphanumeric-character hash based on current timestamp"""
    dt = datetime.now()
    x = int(
        dt.month * 36000 * 24 * 365.25
        + dt.day * 36000 * 24
        + dt.hour * 36000
        + dt.minute * 600
        + dt.second * 10
        + dt.microsecond / 10000
    )
    timehash = ''.join([ALPHANUMERIC_CHARS[v] for v in number_to_base(x, 36)])
    return f'{timehash:>07}'


class SanitisationError(Exception):
    """Something went wrong in Molecule/SMILES sanitisation"""

    ...


class MissingTagError(ValueError):
    """A referenced tag is not in the tag vocabulary.

    Subclasses :class:`ValueError` to match the ``MethodService`` convention for
    "referenced thing is not registered".
    """

    ...


# SQLSTATE for "permission denied" (insufficient_privilege)
INSUFFICIENT_PRIVILEGE = '42501'


class InsufficientPrivilegeError(PermissionError):
    """The database login lacks the rights for a write, typically adding to a
    vocabulary (tags, methods) that only DesignDB admins may extend.

    Rights come from the user's infrastructure login, not from HIPPO; HIPPO only
    attempts the write and explains the refusal.
    """

    ...


@contextmanager
def admin_write(action: str) -> Iterator[None]:
    """Turn a database permission refusal into :class:`InsufficientPrivilegeError`.

    Wrap writes that need DesignDB admin rights. The write is attempted as
    normal; if Postgres refuses it for lack of privilege, the user is told what
    was refused and who to ask. Every other error passes through unchanged.

    :param action: what was being attempted, completing "<action> needs DesignDB
        admin rights", e.g. ``"Registering compound tag 'foo'"``
    :raises InsufficientPrivilegeError: if the database refused the write for
        lack of privilege
    """
    try:
        yield
    except ProgrammingError as exc:
        if getattr(exc.__cause__, 'sqlstate', None) != INSUFFICIENT_PRIVILEGE:
            raise
        user = connection.settings_dict.get('USER') or 'unknown'
        raise InsufficientPrivilegeError(
            f'{action} needs DesignDB admin rights, which your database login '
            f'({user!r}) does not have. Ask a DesignDB administrator to do it '
            'for you.'
        ) from exc


def make_warn_once_per_key():
    """Warn once per field type in sdf file.

    When attribute is defined but broken in all molecules, no need to
    complain every time.

    Instatiate at the beginning of the loading process and pass where
    needed.

    """
    warned = set()

    def warn(key, msg):
        if key not in warned:
            print(f'WARNING: {msg}')
            warned.add(key)

    return warn


# TODO: move
class ScoreSubquery(Subquery):
    def __init__(self, scoring_method):
        # avoiding circular imports
        from .models import ScoreValueModel

        query = ScoreValueModel.objects.filter(
            pose=OuterRef('pk'),
            compound=OuterRef('compound'),
            scoring_method__method_name=scoring_method,
        ).values('score')[:1]
        super().__init__(query)


# Don't understand the distinct here. Shouldn't have to use it.
# Workaround for missing ArrayAgg in sqlite, can get rid of when
# moving to postgres
class JsonGroupArray(Aggregate):
    function = 'json_group_array'
    # template = "%(function)s(%(expressions)s)"
    template = '%(function)s(DISTINCT %(expressions)s)'


def normalize_string_list(x):
    """Convert string representation of list to proper list"""
    if not x:
        return []
    if isinstance(x, list):
        # return list(set(x))
        return x
    if isinstance(x, str):
        # try JSON first
        try:
            parsed = json.loads(x)
            if isinstance(parsed, list):
                # return list(set(parsed))
                return parsed
        except Exception:
            pass

        # fallback for python-style strings
        try:
            parsed = ast.literal_eval(x)
            if isinstance(parsed, list):
                # return list(set(parsed))
                return parsed
        except Exception:
            pass

        # ultimate fallback, comma-separated string
        try:
            splits = x.split(',')
            if isinstance(splits, list):
                return splits
        except Exception:
            pass
    return []


def superparent(mol: Chem.Mol) -> Chem.Mol:
    return rdMolStandardize.SuperParent(mol)


def registration_hash_tautomer_insensitive(mol: Chem.Mol) -> str:
    layers = RegistrationHash.GetMolLayers(
        mol,
        escape='',
        enable_tautomer_hash_v2=True,
    )
    return RegistrationHash.GetMolHash(
        layers,
        RegistrationHash.HashScheme.TAUTOMER_INSENSITIVE_LAYERS,
    )


# Bound on the SMILES -> hash memo. A single load of <=100k molecules never fills
# it; past that we stop inserting rather than evicting, keeping the entries from
# earliest in the run (which is where repeats cluster).
_HASH_CACHE_MAX = 100_000
_HASH_CACHE: dict[str, str] = {}

# Below this many uncached molecules, process-pool startup costs more than the
# parallelism saves, so hashing runs serially.
PARALLEL_HASH_MIN_MOLECULES = 100


# PostgreSQL's wire protocol encodes the bind-parameter count as int16, capping a
# single statement at 65535 parameters. Django only chunks bulk_create for SQLite
# (BaseDatabaseOperations.bulk_batch_size returns len(objs)), so on PostgreSQL an
# unbounded bulk_create of N rows sends N x fields parameters in one statement and
# fails once that exceeds the cap -- around 4,400 rows for PoseModel.
PG_MAX_BIND_PARAMS = 65535

# Independently of the protocol cap, keep single statements to a sane size: pose
# rows carry mol blobs, and a multi-megabyte INSERT is bad for memory and for lock
# duration even when it is legal.
MAX_ROWS_PER_STATEMENT = 5000


# Target wire size for one bulk statement. Both psycopg (building the parameter
# array) and PostgreSQL (parsing the multi-VALUES statement) hold the whole thing
# in memory, so this bounds peak RAM on both ends and keeps lock duration and
# retry cost sane. The latency-vs-size curve is flat across a wide middle range,
# so the exact value matters far less than avoiding the extremes.
TARGET_STATEMENT_BYTES = 4_000_000

# How many rows to sample when estimating row width. Cheap, and spread across the
# batch rather than taken from the front, since SDFs are often size-ordered.
ROW_SIZE_SAMPLE = 64


def _value_bytes(value) -> int:
    """Approximate wire size of a single bound parameter."""
    if value is None:
        return 1
    if isinstance(value, memoryview | bytes | bytearray):
        return len(value)
    if isinstance(value, str):
        return len(value.encode('utf-8', 'ignore'))
    # numbers, booleans, dates: bound as fixed-width values
    return 8


def estimate_row_bytes(model, objs, *, fields: list[str] | None = None) -> float:
    """Estimate the average wire size of one row of ``objs``.

    Values are measured through each field's ``get_prep_value``, so this reflects
    what is actually bound -- notably ``MolField``, which sends
    ``mol.ToBinary()`` rather than the molblock text.

    :param model: the Django model being written
    :param objs: the objects about to be written
    :param fields: field names to measure, or ``None`` for every concrete field.
        Pass the subset for ``bulk_update``, which only binds the updated columns.
    :returns: mean bytes per row, or ``0.0`` if it cannot be estimated
    """
    if not objs:
        return 0.0

    if fields is None:
        model_fields = list(model._meta.concrete_fields)
    else:
        model_fields = [model._meta.get_field(name) for name in fields]

    stride = max(1, len(objs) // ROW_SIZE_SAMPLE)
    sample = objs[::stride][:ROW_SIZE_SAMPLE]
    if not sample:
        return 0.0

    total = 0
    for obj in sample:
        for field in model_fields:
            try:
                value = field.get_prep_value(getattr(obj, field.attname, None))
            except Exception:
                # a field whose prep raises on an unsaved object: fall back to the
                # raw attribute rather than failing the size estimate
                value = getattr(obj, field.attname, None)
            total += _value_bytes(value)

    return total / len(sample)


def safe_batch_size(
    model,
    *,
    objs=None,
    fields: list[str] | None = None,
    requested: int | None = None,
    target_bytes: int = TARGET_STATEMENT_BYTES,
    extra_params: int = 0,
) -> int:
    """Rows to send per ``bulk_create``/``bulk_update`` statement for ``model``.

    Two independent limits apply, and the smaller wins:

    1. **Bind parameters** -- PostgreSQL encodes the parameter count as int16, so
       ``rows x fields`` must stay under 65535. This is a correctness ceiling, not
       a tuning knob.
    2. **Statement bytes** -- estimated from ``objs`` and held near
       ``target_bytes``. This is what actually differs between models: a pose row
       carries a serialised molecule while a junction row is a few integers, so
       equal row counts produce statements orders of magnitude apart.

    Sampling ``objs`` rather than assuming a fixed row width means unusually large
    ligands automatically produce smaller batches instead of oversized statements.

    Without ``objs`` only the parameter ceiling applies (bounded by
    :data:`MAX_ROWS_PER_STATEMENT`), which is safe but ignores row width.

    :param model: the Django model being written
    :param objs: the objects about to be written, sampled to estimate row width
    :param fields: field subset being written, for ``bulk_update``
    :param requested: caller's explicit preference; honoured up to the parameter
        ceiling and clamped above it
    :param target_bytes: desired wire size per statement
    :param extra_params: additional per-statement parameters to leave room for
    :returns: rows per statement, at least 1
    """
    n_fields = max(1, len(fields or model._meta.concrete_fields))
    # 5% headroom: the field count is an estimate of what Django actually binds,
    # and landing exactly on the cap leaves no room for a conflict clause or an
    # extra expression to tip it over.
    budget = int((PG_MAX_BIND_PARAMS - extra_params) * 0.95)
    param_ceiling = max(1, budget // n_fields)

    if requested is not None:
        if requested > param_ceiling:
            mrich.warning(
                f'batch_size={requested} exceeds the safe maximum for '
                f'{model.__name__} ({param_ceiling} rows x {n_fields} fields vs '
                f'the {PG_MAX_BIND_PARAMS}-parameter cap); using {param_ceiling}'
            )
            return param_ceiling
        return max(1, requested)

    ceiling = min(param_ceiling, MAX_ROWS_PER_STATEMENT)

    row_bytes = estimate_row_bytes(model, objs, fields=fields) if objs else 0.0
    if row_bytes <= 0:
        return ceiling

    return max(1, min(int(target_bytes / row_bytes), ceiling))


def default_hash_workers() -> int:
    """Default worker count for parallel registration hashing."""
    return min(8, os.cpu_count() or 1)


def _hash_worker(smiles: str) -> tuple[str, str | None]:
    """Compute one registration hash. Runs in a *separate process*.

    .. warning::
       Must stay pure RDKit. This is forked from a parent that may hold an open
       database connection inside a transaction; forking is safe only as long as
       the child never touches that socket. Do not import or call anything
       Django/ORM here.

    :param smiles: SMILES to hash
    :returns: ``(smiles, hash)``, or ``(smiles, None)`` if it could not be hashed
    """
    try:
        return smiles, _compute_compound_hash(smiles)
    except ValueError:
        return smiles, None


def _compute_compound_hash(smiles: str) -> str:
    """Uncached, unparallelised registration hash for a SMILES.

    :param smiles: SMILES to hash
    :returns: tautomer-insensitive registration hash
    :raises ValueError: if the SMILES cannot be parsed or SuperParent fails
    """
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        raise ValueError(f'Could not parse SMILES: {smiles!r}')
    try:
        sp = rdMolStandardize.SuperParent(mol)
    except Exception as e:
        raise ValueError(f'SuperParent failed: {e}') from e
    return registration_hash_tautomer_insensitive(sp)


def compound_hash_from_smiles(smiles: str) -> str:
    """Registration hash for a SMILES, via its SuperParent.

    ``rdMolStandardize.SuperParent`` is the dominant pure-CPU cost of SDF
    ingestion (~6 ms/mol), so results are memoised on the SMILES string. The memo
    only helps when a SMILES repeats; for a batch of mostly-distinct molecules use
    :func:`compound_hashes_from_smiles`, which parallelises instead.

    :param smiles: SMILES to hash
    :returns: tautomer-insensitive registration hash
    :raises ValueError: if the SMILES cannot be parsed or SuperParent fails
    """
    cached = _HASH_CACHE.get(smiles)
    if cached is not None:
        return cached
    result = _compute_compound_hash(smiles)
    if len(_HASH_CACHE) < _HASH_CACHE_MAX:
        _HASH_CACHE[smiles] = result
    return result


def compound_hashes_from_smiles(
    smiles_list: list[str],
    *,
    max_workers: int | None = None,
) -> dict[str, str]:
    """Registration hashes for many SMILES, in parallel.

    SuperParent is CPU-bound C++ that does not release the GIL, so threads cannot
    help; separate processes scale close to linearly (~6x on 8 workers for a
    1000-molecule SDF). Only cache misses are dispatched, and only SMILES strings
    cross the process boundary -- never RDKit molecules or ORM objects.

    Falls back to serial hashing when there is little to do or when a pool cannot
    be started, so this is always safe to call.

    :param smiles_list: SMILES to hash; duplicates are collapsed
    :param max_workers: worker processes, defaulting to
        :func:`default_hash_workers` (``min(8, cpu_count())``). ``1`` forces serial.
    :returns: mapping of SMILES to hash. SMILES that could not be hashed are
        absent, matching :func:`compound_hash_from_smiles` raising for them.
    """
    result: dict[str, str] = {}
    todo: list[str] = []

    for smiles in dict.fromkeys(smiles_list):
        cached = _HASH_CACHE.get(smiles)
        if cached is not None:
            result[smiles] = cached
        else:
            todo.append(smiles)

    if not todo:
        return result

    workers = default_hash_workers() if max_workers is None else max_workers

    if workers > 1 and len(todo) >= PARALLEL_HASH_MIN_MOLECULES:
        # ~4 chunks per worker: enough to even out molecules of differing cost
        # without paying per-task IPC on every molecule
        chunksize = max(1, len(todo) // (workers * 4))
        try:
            with ProcessPoolExecutor(max_workers=workers) as pool:
                pairs = list(pool.map(_hash_worker, todo, chunksize=chunksize))
        except Exception as e:
            # restricted//sandboxed environments may refuse to fork
            mrich.warning(f'Parallel hashing unavailable ({e}); falling back')
            pairs = [_hash_worker(smiles) for smiles in todo]
    else:
        pairs = [_hash_worker(smiles) for smiles in todo]

    for smiles, value in pairs:
        if value is None:
            continue
        result[smiles] = value
        if len(_HASH_CACHE) < _HASH_CACHE_MAX:
            _HASH_CACHE[smiles] = value

    return result


def profile(output_file='profile.prof'):
    """Function profiler decorator.

    Usage: just add the decorator
    @profile(<filename>)
    """

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):

            profiler = cProfile.Profile()
            profiler.enable()
            try:
                return func(*args, **kwargs)
            finally:
                profiler.disable()
                profiler.dump_stats(output_file)

        return wrapper

    return decorator
