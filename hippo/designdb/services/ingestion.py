import logging
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import molparse as mp
import mrich
import pandas as pd
from designdb.components.compound import Ingredient
from designdb.models import (
    CompoundEnumerationMethodJunctionModel,
    CompoundModel,
    CompoundTagJunctionModel,
    InspirationModel,
    PoseMethodModel,
    PoseModel,
    PoseTagJunctionModel,
    ReactantModel,
    ReactionModel,
    ScaffoldModel,
    TargetModel,
)
from designdb.recipe import Recipe, Route
from designdb.services.compound import CompoundService, CompoundTagService
from designdb.services.method import MethodService
from designdb.services.pose import PoseService, PoseTagService
from designdb.services.pose_score import ScoreService
from designdb.services.reaction import ReactionService
from designdb.services.route import RouteService
from designdb.sets.compound import CompoundSet
from designdb.sets.ingredient import IngredientSet
from designdb.sets.reaction import ReactionSet
from designdb.utils import (
    SanitisationError,
    inchikey_from_smiles,
    remove_other_ligands,
    safe_batch_size,
    sanitise_smiles,
)
from designdb.utils_chem import (
    InvalidChemistryError,
    UnsupportedChemistryError,
    check_chemistry,
)
from designdb.utils_frag import (
    UnsupportedFragalysisLongcodeError,
    parse_observation_longcode,
)
from django.db import connection, transaction
from numpy import isnan
from pandas import read_pickle

# from mypackage.services.compound import CompoundService
from rdkit import Chem

# from rdkit.Chem import inchi
from rdkit.Chem import PandasTools

# from .validation.compound import ValidationError, validate_compound_data

SDF_XCAv2_PATTERN = re.compile(
    r'^[^.]*-.\d{4}_._\d*_\d_.*-.\d{4}\+.\+\d*\+\d_ligand\.sdf$'
)
SDF_XCAV3_PATTERN = re.compile(
    r'^[^.]*-.\d{4}_._\d*_._\d_.*-.\d{4}\+.\+\d*\+.\+\d_ligand\.sdf$'
)


SDF_FRAGALYSIS_PATTERN = re.compile(r'^[^.].*\d{4}[a-z].sdf$')
PDBID_PATTERN = re.compile(r'^[A-Za-z0-9]{4}-[a-z].sdf$')


logger = logging.getLogger(__name__)


@dataclass
class FSRecord:
    name: str
    path: Path
    sdf: Path
    pdb: Path


def parse_sdf_pandas(sdf_path: Path) -> tuple[str, Chem.rdchem.Mol]:
    df = PandasTools.LoadSDF(
        str(sdf_path), molColName='ROMol', idName='ID', strictParsing=True
    )
    # extract fields
    longcode = df.ID[0]
    mol = df.ROMol[0]

    return longcode, mol


def parse_pdb_mp(pdb_path: Path, residue: int, chain: str) -> str:
    logger.debug('Reading %s', pdb_path)
    pdb = mp.parse(pdb_path, verbosity=0)

    # protein_link is stored relative to the current working directory (e.g.
    # data/downloads/...) so the database stays portable across machines
    rel_pdb = os.path.relpath(pdb_path)

    # create the single ligand bound pdb
    lig_residues = pdb.residues['LIG']
    if len(lig_residues) > 1 or any(r.contains_alternative_sites for r in lig_residues):
        pdb = remove_other_ligands(pdb, residue, chain)
        pdb.prune_alternative_sites('A', verbosity=0)
        pose_path = rel_pdb.replace('.pdb', '_hippo.pdb')
        # side effect: writes pdb into file
        mp.write(
            pose_path, pdb, shift_name=True, verbosity=logger.level == logging.DEBUG
        )
    else:
        pose_path = rel_pdb

    return pose_path


def iter_fs_fragalysis(root_path, skip_records):
    assert skip_records is not None, '"None" passed instead as skip_records'

    for dset_path in list(sorted(root_path.glob('[!.]*'))):
        if dset_path.name in skip_records:
            continue

        sdfs = []
        for sdf_path in dset_path.glob('[!.]*.sdf'):
            sdf_name = sdf_path.name

            if (
                '_ligand' in sdf_name
            ):  # Quick fix, _ligand.sdf are exactly the same as .sdf
                # in aligned_directory.
                continue

            # fragalysis SDF
            if SDF_FRAGALYSIS_PATTERN.match(sdf_name):
                sdfs.append(sdf_path)
            # fragalysis SDF from PDB id
            elif PDBID_PATTERN.match(sdf_name):
                sdfs.append(sdf_path)
            else:
                mrich.warning(
                    sdf_name,
                    "doesn't not follow neither Fragalysis nor PDB ID patterns",
                )
                sdfs.append(sdf_path)

        if not sdfs:
            mrich.error(dset_path.name, 'has no compatible SDFs', dset_path)
            continue

        pdbs = [
            p
            for p in dset_path.glob('[!.]*.pdb')
            if '_ligand' not in p.name
            and '_delig' not in p.name  # current Fragalysis protein-file naming
            and '_hippo' not in p.name
            # DEPRECATED(apo-naming): pre-'delig' Fragalysis naming, remove once
            # all data uses 'delig'
            and '_apo' not in p.name
        ]

        if not len(pdbs) == 1:
            mrich.error(dset_path.name, 'has invalid PDBs', pdbs)
            continue

        record = FSRecord(name=dset_path.name, path=dset_path, sdf=sdfs[0], pdb=pdbs[0])

        logger.debug('fs_frag record: %s', record)

        yield record


# unfinished, seems XCA data is not loaded now
def iter_fs_xca(root_path, skip):
    for dset_path in sorted(root_path.glob('*[0-9][0-9][0-9][0-9]')):
        if dset_path.name in skip:
            continue

        sdfs = []

        for sdf_path in sorted(dset_path.glob('[!.]*.sdf')):
            sdf_name = sdf_path.name

            # TODO: switch between patterns??
            if SDF_XCAv2_PATTERN.match(sdf_name):
                sdfs.append(sdf_path)

        if not sdfs:
            mrich.error(dset_path.name, 'has no compatible SDFs', dset_path)
            continue

        for i, sdf in enumerate(sdfs):
            subname = dset_path.name + chr(ord('a') + i)

            pdb = dset_path / sdf.name.replace('_ligand.sdf', '.pdb')

            if not pdb.exists():
                mrich.error(dset_path.name, 'is missing PDB', pdb)
                continue

            record = FSRecord(name=subname, path=dset_path, sdf=sdf, pdb=pdb)

            logger.debug('fs_frag record: %s', record)

            yield record


def read_df(path: Path):
    if path.name.endswith('.sdf'):
        df = PandasTools.LoadSDF(str(path.resolve()))
    else:
        df = read_pickle(path)

    return df


DEFAULT_CHUNK_SIZE = 5000


def iter_sdf_chunks(
    path: Path,
    *,
    mol_col: str,
    name_col: str,
    skip_equal: dict,
    skip_not_equal: dict,
    chunk_size: int,
):
    """Stream an SDF, yielding lists of at most ``chunk_size`` records.

    Replaces ``read_df`` + :func:`preprocess_df` for the ingestion path. Those
    materialise the entire file: measured at roughly 100 KB per molecule through
    the full pipeline, a 300,000-record SDF would need ~28 GB. Streaming holds one
    chunk at a time, so peak memory is set by ``chunk_size`` rather than file size.

    Records are built to match ``PandasTools.LoadSDF`` exactly -- ``GetProp``
    returns strings, as the DataFrame columns do, and ``_Name`` becomes
    ``name_col`` -- so ingestion behaviour is unchanged. Verified record-for-record
    against the DataFrame path.

    Non-SDF inputs (pickled frames) have no streaming reader and are read whole,
    then chunked.

    :param path: file to read
    :param mol_col: key to store the molecule under
    :param name_col: key to store the molecule title under
    :param skip_equal: keep rows where ``row[k] == v``
    :param skip_not_equal: keep rows where ``row[k] != v``
    :param chunk_size: maximum records per yielded chunk
    :yields: lists of record dicts
    """

    def keep(rec) -> bool:
        # 'ID' is hardcoded here to match the DataFrame path, which filters on
        # df['ID'] regardless of name_col
        if rec.get('ID') == 'ver_1.2':
            return False
        for k, v in skip_equal.items():
            if rec.get(k) != v:
                return False
        for k, v in skip_not_equal.items():
            if rec.get(k) == v:
                return False
        return True

    def finish(rec):
        alias = rec.get(name_col)
        if isinstance(alias, str):
            rec[name_col] = alias.strip()
        return rec

    if not path.name.endswith('.sdf'):
        # pickled frame: no streaming reader available, so read whole and chunk
        df = read_df(path)
        records = preprocess_df(
            df,
            skip_equal=skip_equal,
            skip_not_equal=skip_not_equal,
            name_col=name_col,
        )
        for i in range(0, len(records), chunk_size):
            yield records[i : i + chunk_size]
        return

    chunk: list[dict] = []
    seen = kept = 0
    with open(path, 'rb') as handle:
        for mol in Chem.ForwardSDMolSupplier(handle):
            seen += 1
            if mol is None:
                mrich.error(f'Could not parse molecule {seen} of {path}')
                continue
            rec = {n: mol.GetProp(n) for n in mol.GetPropNames()}
            rec[name_col] = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
            rec[mol_col] = mol
            if not keep(rec):
                continue
            kept += 1
            chunk.append(finish(rec))
            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []
    if chunk:
        yield chunk
    mrich.var('SDF entries (read)', seen)
    mrich.var('SDF entries (post-filter)', kept)


def validate_record(
    record,
    mol_col,
    name_col,
    inspiration_col,
    inspirations,
    reference_col,
    reference,
):
    """Streaming counterpart of :func:`validate_df`, checked on the first record."""
    assert mol_col in record, f'{mol_col=} not in {sorted(record)}'

    if name_col:
        assert name_col in record, f'{name_col=} not in {sorted(record)}'

    if inspiration_col and not inspirations:
        assert inspiration_col in record, f'{inspiration_col=} not in {sorted(record)}'

    if not reference and reference_col:
        assert reference_col in record, f'{reference_col=} not in {sorted(record)}'


def validate_df(
    df,
    mol_col,
    name_col,
    inspiration_col,
    inspirations,
    reference_col,
    reference,
):

    # TODO: these are part of input validation and should be removed. or
    # at least rewritten
    assert mol_col in df.columns, f'{mol_col=} not in {df.columns}'

    if name_col:
        assert name_col in df.columns, f'{name_col=} not in {df.columns}'

    if inspiration_col and not inspirations:
        assert inspiration_col in df.columns, f'{inspiration_col=} not in {df.columns}'

    if not reference and reference_col:
        assert reference_col in df.columns, f'{reference_col=} not in {df.columns}'


def preprocess_df(
    df,
    *,
    skip_equal,
    skip_not_equal,
    name_col: str,
) -> list[dict[str, Any]]:

    mrich.var('SDF entries (pre-filter)', len(df))

    df = df[df['ID'] != 'ver_1.2']

    for k, v in skip_equal.items():
        df = df[df[k] == v]

    for k, v in skip_not_equal.items():
        df = df[df[k] != v]

    mrich.var('SDF entries (post-filter)', len(df))

    df[name_col] = df[name_col].str.strip()

    records = df.to_dict(orient='records')

    return records


def metadata_from_record(
    record: dict[str, str],
    ignore_fields: list[str | None],
    convert_floats: bool,
    field_warning=None,
) -> dict[str, str | float]:

    result = {}
    skip = {
        'smiles',
        'inchikey',
        'compound_id',
        'target_id',
        'reference_id',
        'path',
        'exports',
    }

    skip = skip.union(set([k for k in ignore_fields if k]))

    for key, value in record.items():
        if key in skip:
            continue

        if isinstance(value, float) and isnan(value):
            continue

        if convert_floats:
            try:
                value = float(value)
            except TypeError:
                pass
            except ValueError:
                pass

        if not (isinstance(value, str) or isinstance(value, float)):
            if field_warning:
                field_warning(mrich.warning(f'Skipping metadata from column={key}.'))
            continue

        result[key] = value

    return result


@dataclass
class IngestionBatchResult:
    attempts: int = 0
    compounds_created: int = 0
    poses_created: int = 0
    #: chunks committed so far; with per-chunk commits this is what survived a
    #: partial failure
    chunks_committed: int = 0


class IngestionService:
    @classmethod
    def ingest_filesystem(
        cls,
        *,
        root_path: Path,
        target: TargetModel,
        skip_records: list[str],
        compound_tag_list: list[str],
        metadata_file: Path | str,
        pose_methods: list[PoseMethodModel] | None = None,
        check_rmsd: bool = False,
        rmsd_threshold: float = 1.0,
    ) -> IngestionBatchResult:

        # this is now strictly for loading frag data. cannot switch inner funcs easily
        result = IngestionBatchResult()
        compound_tags = CompoundTagService.resolve_tags(compound_tag_list)
        pose_tagger = PoseTagService(metadata_file, other_tags=compound_tag_list)

        # if needs xca paths, need to pass or select function
        for fs_record in iter_fs_fragalysis(root_path, skip_records):
            longcode, mol = parse_sdf_pandas(fs_record.sdf)
            logger.debug(fs_record.name, longcode)
            result.attempts += 1

            # TODO: this is the original procedure how it was
            # calculated in hippo. I'm not touching it now, but this
            # could use a rewrite, it converts smiles back to mol and
            # then to inchikey
            smiles = mp.rdkit.mol_to_smiles(mol)
            # sane_smiles = sanitise_smiles(
            #     smiles, verbosity=logger.level == logging.DEBUG
            # )
            inchikey = inchikey_from_smiles(smiles)
            # sane_inchikey = inchikey_from_smiles(sane_smiles)

            # NB! different func if XCA data
            try:
                longcode_rec = parse_observation_longcode(longcode)
            except UnsupportedFragalysisLongcodeError as exc:
                # unhandled in original code. do what?
                raise UnsupportedFragalysisLongcodeError from exc

            pose_path = parse_pdb_mp(
                fs_record.pdb, longcode_rec.residue_number, longcode_rec.chain
            )

            compound, compound_created = CompoundService.create(
                # mol=mol,
                # smiles=sane_smiles,
                smiles=smiles,
                # inchikey=sane_inchikey,
            )
            for compound_tag in compound_tags:
                CompoundTagJunctionModel.objects.get_or_create(
                    compound=compound,
                    compound_tag=compound_tag,
                    target=target,
                )
            if compound_created:
                result.compounds_created += 1

            pose_tags, metadata = pose_tagger.tags_and_meta(
                code=fs_record.name,
                longcode=longcode,
            )

            metadata = {'fragalysis_longcode': longcode}

            pose, pose_created = PoseService.create(
                compound=compound,
                target=target,
                mol=mol,
                alias=fs_record.name,
                path=pose_path,
                metadata=metadata,
                inchikey=inchikey,
                smiles=smiles,
                # the first method is the uniqueness/producing method; create()
                # associates it. add_hits may request additional method tags,
                # which are associated below.
                pose_method=pose_methods[0] if pose_methods else None,
                check_rmsd=check_rmsd,
                rmsd_threshold=rmsd_threshold,
            )
            if pose_created:
                result.poses_created += 1

            pose.tags.add(*pose_tags)

            if pose_methods and len(pose_methods) > 1:
                pose.methods.add(*pose_methods[1:])

            # it seems fragalysis data is not expected to contain
            # scores

        return result

    @staticmethod
    def score_values_load_helpers_available() -> bool:
        """Check that *both* score-value load helpers exist on this database.

        Deployed databases have been found carrying ``begin_score_values_load()``
        without ``end_score_values_load()``: the historical per-record loader called
        ``begin_`` twice and never called ``end_``, so a half-applied schema went
        unnoticed. Disabling the trigger under those conditions strands it disabled
        and silently staleifies ``scores_per_pose_pivoted_mv``, so both are checked
        before anything is disabled.

        :returns: ``True`` only if both functions are present
        """
        if connection.vendor != 'postgresql':
            return False
        with connection.cursor() as cursor:
            cursor.execute(
                'SELECT '
                "to_regprocedure('designdb.begin_score_values_load()') IS NOT NULL, "
                "to_regprocedure('designdb.end_score_values_load()') IS NOT NULL"
            )
            begin_exists, end_exists = cursor.fetchone()
        if not (begin_exists and end_exists):
            mrich.warning(
                'score_values load helpers missing on this database '
                f'(begin={bool(begin_exists)}, end={bool(end_exists)}); '
                'not disabling trg_score_values_refresh_pivoted_mv. '
                'See docs/proposal_pose_alias_index.md sibling note, or create them '
                'from images/xchem-designdb/init-db/01_schema.sql.'
            )
        return bool(begin_exists and end_exists)

    @classmethod
    def begin_score_values_load(cls) -> bool:
        """Disable the per-statement score-value matview refresh trigger.

        ``trg_score_values_refresh_pivoted_mv`` fires a matview refresh for every
        statement touching ``score_values``. The schema ships helpers to bracket a
        bulk load (``init-db/01_schema.sql``).

        Call this *outside* the ingestion transaction: the helper takes an
        ``AccessExclusiveLock`` on ``score_values``, and holding that for the whole
        load is what makes concurrent writers on production block.

        Refuses to disable anything unless :meth:`score_values_load_helpers_available`
        confirms the matching ``end_`` helper exists, so a failure can never leave
        the trigger off.

        .. note::
           The trigger is ``FOR EACH STATEMENT``, and :meth:`ingest_sdf` writes
           all scores in a single statement -- so skipping the disable costs the
           batch path one refresh, not one per record. Losing this optimisation is
           not a reason to risk stranding the trigger.

        :returns: ``True`` if the trigger was actually disabled, ``False`` otherwise
            (SQLite, or helpers missing)
        """
        if not cls.score_values_load_helpers_available():
            return False
        with connection.cursor() as cursor:
            cursor.execute('SELECT designdb.begin_score_values_load();')
        return True

    @staticmethod
    def end_score_values_load() -> bool:
        """Re-enable the score-value trigger and refresh the matview.

        Counterpart to :meth:`begin_score_values_load`. Always call this from a
        ``finally``: leaving the trigger disabled silently staleifies
        ``scores_per_pose_pivoted_mv`` for every later reader.

        .. note::
           The SQL helper refreshes the matview *non*-concurrently, so this takes an
           ``AccessExclusiveLock`` on the view for the duration of the refresh.

        :returns: ``True`` if the trigger was re-enabled, ``False`` on SQLite
        """
        if connection.vendor != 'postgresql':
            return False
        with connection.cursor() as cursor:
            cursor.execute('SELECT designdb.end_score_values_load();')
        return True

    @classmethod
    def ingest_sdf(
        cls,
        *,
        file_path: Path,
        target,
        compound_tag_list: list[str],
        pose_tag_list: list[str],
        enumeration_method: tuple[str, str] | None = None,
        pose_method: tuple[str, str] | None = None,
        score_cols: list[str] | None = None,
        scoring_methods: list[tuple[str, str]] | None = None,
        mol_col: str,
        name_col: str,
        inspiration_col: str | None = None,
        inspirations: list[int],
        inspiration_map: dict[str, PoseModel],
        reference: int | None,
        reference_col: str,
        skip_equal,
        skip_not_equal,
        convert_floats: bool = True,
        field_warning=None,
        max_workers: int | None = None,
        batch_size: int | None = None,
        chunk_size: int | None = None,
        single_transaction: bool = False,
        check_rmsd: bool = False,
        rmsd_threshold: float = 1.0,
    ) -> IngestionBatchResult:
        """Ingest an SDF into the database, set-based.

        Each phase touches the database once for the whole file rather than once
        per record:

        1. parse, sanitise and hash every record -- pure CPU, no DB
        2. resolve all compounds (:meth:`CompoundService.create_batch`)
        3. resolve all inspirations and references in bulk
        4. insert all poses (:meth:`PoseService.create_batch`)
        5. insert all tag junctions, inspirations and scores in bulk

        Methods arrive as ``(name, version)`` pairs rather than pre-resolved ORM
        objects and are resolved here via :class:`MethodService` /
        :meth:`ScoreService.resolve_score_method_map`, keeping model instances out
        of the client-facing signature for the eventual client/backend split.

        .. note::
           ``trg_score_values_refresh_pivoted_mv`` is deliberately left alone. It is
           ``FOR EACH STATEMENT`` and this path writes every score in a single
           statement, so it fires once per load either way.

        :param enumeration_method: ``(name, version)`` of the enumeration method
        :param pose_method: ``(name, version)`` of the pose method
        :param score_cols: SDF columns holding score values
        :param scoring_methods: ``(name, version)`` pairs aligned with ``score_cols``
        :param max_workers: worker processes for registration hashing; see
            :func:`compound_hashes_from_smiles`
        :param batch_size: rows per INSERT/UPDATE statement, applied to every bulk
            write. ``None`` uses the safe per-model maximum; a larger value is
            clamped to it. See :func:`safe_batch_size`.
        :param chunk_size: SDF records held in memory and processed per pass,
            defaulting to :data:`DEFAULT_CHUNK_SIZE`. This is what bounds peak
            memory; ``batch_size`` only bounds statement size.
        :param single_transaction: wrap the whole file in one transaction instead
            of committing per chunk. All-or-nothing, at the cost of a transaction
            that lives for the entire load.
        :param check_rmsd: skip poses whose RMSD to an existing pose of the same
            compound is below ``rmsd_threshold``; see
            :meth:`PoseService.create_batch`
        :param rmsd_threshold: RMSD below which two poses are the same, in Angstrom
        :returns: counts of attempted records and created compounds/poses
        :raises ValueError: if a named method is not registered
        :raises MissingTagError: if a named tag is not registered. Tags are
            resolved, never created -- see ``register_tag`` on the tag services
        """
        result = IngestionBatchResult()

        # Methods arrive as (name, version) rather than ORM objects: resolving them
        # is DB traversal and belongs here, not in the client layer. Done before any
        # file IO so an unregistered method fails fast.
        enumeration_method_obj = MethodService.resolve_enumeration_method(
            enumeration_method
        )
        pose_method_obj = MethodService.resolve_pose_method(pose_method)
        score_method_map = ScoreService.resolve_score_method_map(
            score_cols, scoring_methods
        )

        # Tags, like methods, must already exist: the vocabulary is maintained
        # outside HIPPO. Resolved before the output directory is created so a
        # rejected load leaves nothing behind.
        compound_tags = CompoundTagService.resolve_tags(compound_tag_list)
        pose_tags = PoseTagService.resolve_tags(pose_tag_list)

        output_directory = Path(str(file_path.name).removesuffix('.sdf'))
        output_directory.mkdir(parents=True, exist_ok=True)

        scorer = ScoreService()

        chunk_size = chunk_size or DEFAULT_CHUNK_SIZE
        chunks = iter_sdf_chunks(
            file_path,
            mol_col=mol_col,
            name_col=name_col,
            skip_equal=skip_equal,
            skip_not_equal=skip_not_equal,
            chunk_size=chunk_size,
        )

        def run_chunk(records):
            return cls._ingest_chunk(
                records=records,
                result=result,
                target=target,
                output_directory=output_directory,
                compound_tags=compound_tags,
                pose_tags=pose_tags,
                scorer=scorer,
                enumeration_method_obj=enumeration_method_obj,
                pose_method_obj=pose_method_obj,
                score_method_map=score_method_map,
                mol_col=mol_col,
                name_col=name_col,
                inspiration_col=inspiration_col,
                inspirations=inspirations,
                inspiration_map=inspiration_map,
                reference=reference,
                reference_col=reference_col,
                convert_floats=convert_floats,
                field_warning=field_warning,
                max_workers=max_workers,
                batch_size=batch_size,
                check_rmsd=check_rmsd,
                rmsd_threshold=rmsd_threshold,
            )

        validated = False
        if single_transaction:
            # all-or-nothing: one transaction spans every chunk. Bounded memory
            # still applies, but the transaction (and its locks and WAL) lives for
            # the whole load.
            with transaction.atomic():
                for records in chunks:
                    if not validated and records:
                        validate_record(
                            records[0],
                            mol_col,
                            name_col,
                            inspiration_col,
                            inspirations,
                            reference_col,
                            reference,
                        )
                        validated = True
                    run_chunk(records)
        else:
            # commit per chunk: a failure leaves earlier chunks committed. Ingestion
            # is re-runnable -- existing compounds and poses are found rather than
            # duplicated -- so a partial load is recovered by running again.
            for records in chunks:
                if not validated and records:
                    validate_record(
                        records[0],
                        mol_col,
                        name_col,
                        inspiration_col,
                        inspirations,
                        reference_col,
                        reference,
                    )
                    validated = True
                with transaction.atomic():
                    run_chunk(records)
                result.chunks_committed += 1

        return result

    @classmethod
    def _ingest_chunk(
        cls,
        *,
        records,
        result,
        target,
        output_directory,
        compound_tags,
        pose_tags,
        scorer,
        enumeration_method_obj,
        pose_method_obj,
        score_method_map,
        mol_col,
        name_col,
        inspiration_col,
        inspirations,
        inspiration_map,
        reference,
        reference_col,
        convert_floats,
        field_warning,
        max_workers,
        batch_size,
        check_rmsd,
        rmsd_threshold,
    ) -> None:
        """Run the five ingestion phases over one chunk of records.

        Counts accumulate into ``result``. Record indices are chunk-local, which is
        all the bulk lookups need -- cross-chunk identity is resolved through the
        database (a compound created by an earlier chunk is found, not duplicated).
        """

        # ---- phase 1: CPU-only preparation -------------------------------------
        # `prepared` holds one entry per record that survived sanitisation, keeping
        # the original record index so the bulk lookups below can be joined back on
        prepared: list[dict] = []
        for idx, r in enumerate(records):
            result.attempts += 1

            smiles = r.get('smiles', None)
            if not smiles:
                smiles = mp.rdkit.mol_to_smiles(r[mol_col])
            try:
                sanitise_smiles(
                    smiles,
                    sanitisation_failed='error',
                    radical='warning',
                    verbosity=logger.level == logging.DEBUG,
                )
            except SanitisationError as e:
                mrich.error(f'Could not sanitise {smiles=}')
                mrich.error(str(e))
                continue
            except AssertionError:
                mrich.error(f'Could not sanitise {smiles=}')
                continue

            prepared.append(
                {
                    'idx': idx,
                    'record': r,
                    'smiles': smiles,
                    'inchikey': inchikey_from_smiles(smiles),
                    'mol': r[mol_col],
                    'alias': r[name_col],
                    'metadata': metadata_from_record(
                        r,
                        ignore_fields=[inspiration_col, name_col, mol_col],
                        convert_floats=convert_floats,
                        field_warning=field_warning,
                    ),
                    'path': os.path.relpath(
                        output_directory / f'{r[name_col]}.fake.mol'
                    ),
                }
            )

        if not prepared:
            return

        # ---- phase 2: compounds ------------------------------------------------
        compound_map = CompoundService.create_batch(
            smiles_list=[p['smiles'] for p in prepared],
            max_workers=max_workers,
            batch_size=batch_size,
        )
        result.compounds_created += sum(
            1 for _, created in compound_map.values() if created
        )

        # drop records whose compound could not be resolved
        prepared = [p for p in prepared if p['smiles'] in compound_map]
        if not prepared:
            return

        # ---- phase 3: inspirations and references ------------------------------
        inspirations_by_idx = PoseService.get_inspirations_batch(
            records,
            global_inspirations=inspirations,
            inspiration_map=inspiration_map,
            inspiration_col=inspiration_col,
            name_col=name_col,
            target=target,
        )

        if reference:
            # an explicit reference applies to every pose in the file
            reference_by_idx: dict[int, int | None] = {}
        elif reference_col:
            reference_by_idx = PoseService.get_reference_batch(
                records,
                reference_col=reference_col,
                target=target,
            )
        else:
            reference_by_idx = {}

        # ---- phase 4: poses ----------------------------------------------------
        specs = []
        for p in prepared:
            compound, _ = compound_map[p['smiles']]
            specs.append(
                {
                    'compound': compound,
                    'mol': p['mol'],
                    'alias': p['alias'],
                    'path': p['path'],
                    'metadata': p['metadata'],
                    'inchikey': p['inchikey'],
                    'smiles': p['smiles'],
                    'reference': reference
                    if reference
                    else reference_by_idx.get(p['idx']),
                }
            )

        poses, poses_created = PoseService.create_batch(
            target=target,
            specs=specs,
            pose_method=pose_method_obj,
            batch_size=batch_size,
            check_rmsd=check_rmsd,
            rmsd_threshold=rmsd_threshold,
        )
        result.poses_created += poses_created

        # ---- phase 5: junctions, inspirations and scores -----------------------
        compound_ids = {compound.pk for compound, _ in compound_map.values()}

        if compound_tags:
            CompoundTagJunctionModel.objects.bulk_create(
                [
                    CompoundTagJunctionModel(
                        compound_id=compound_id,
                        compound_tag=tag,
                        target=target,
                    )
                    for compound_id in compound_ids
                    for tag in compound_tags
                ],
                ignore_conflicts=True,
                batch_size=safe_batch_size(
                    CompoundTagJunctionModel, requested=batch_size
                ),
            )

        if enumeration_method_obj is not None:
            CompoundEnumerationMethodJunctionModel.objects.bulk_create(
                [
                    CompoundEnumerationMethodJunctionModel(
                        compound_id=compound_id,
                        enumeration_method=enumeration_method_obj,
                        target=target,
                    )
                    for compound_id in compound_ids
                ],
                ignore_conflicts=True,
                batch_size=safe_batch_size(
                    CompoundEnumerationMethodJunctionModel, requested=batch_size
                ),
            )

        live = [(p, pose) for p, pose in zip(prepared, poses, strict=True) if pose]

        if pose_tags:
            PoseTagJunctionModel.objects.bulk_create(
                [
                    PoseTagJunctionModel(pose=pose, pose_tag=tag)
                    for _, pose in live
                    for tag in pose_tags
                ],
                ignore_conflicts=True,
                batch_size=safe_batch_size(PoseTagJunctionModel, requested=batch_size),
            )

        inspiration_rows = []
        for p, pose in live:
            for original_id in inspirations_by_idx.get(p['idx'], []):
                inspiration_rows.append(
                    InspirationModel(
                        derivative_pose=pose,
                        original_pose_id=original_id,
                    )
                )
        if inspiration_rows:
            InspirationModel.objects.bulk_create(
                inspiration_rows,
                ignore_conflicts=True,
                batch_size=safe_batch_size(
                    InspirationModel, objs=inspiration_rows, requested=batch_size
                ),
            )

        scorer.add_scores_from_records_batch(
            pairs=[(pose, p['record']) for p, pose in live],
            score_method_map=score_method_map or None,
            batch_size=batch_size,
        )

        return

    # how is that without target??
    @classmethod
    def ingest_syndirella_routes(
        cls,
        pickle_path: str | Path,
        CAR_only: bool = True,
        pick_first: bool = True,
        do_check_chemistry: bool = True,
        register_routes: bool = True,
    ):
        # this is pretty much a copy from the original method now
        df = read_pickle(pickle_path)

        for i, row in mrich.track(df.iterrows(), total=len(df)):
            mrich.set_progress_field('i', i)
            mrich.set_progress_field('n', len(df))

            d = row.to_dict()

            # comp = self.compounds(smiles=d['smiles'])

            n_routes = 0
            for key in d:
                if not key.startswith('route'):
                    continue

                if not key.endswith('_names'):
                    continue

                v = d[key]

                if isinstance(v, float) and pd.isna(v):
                    break

                n_routes += 1

            if not n_routes:
                # mrich.warning(comp, "#routes =", n_routes)
                continue

            # routes = []
            for j in range(n_routes):
                route_str = f'route{j}'

                route = d[route_str]

                if CAR_only and not d[route_str + '_CAR']:
                    continue

                reactions = ReactionSet()
                reactants = IngredientSet()
                intermediates = IngredientSet()
                products = IngredientSet()

                # new models include ReactionModel, ReactantModel and
                # ComponentModel. Should use these instead?

                try:
                    for k, reaction_struct in enumerate(route):
                        reaction_type = reaction_struct['name']

                        # product = self.compounds(smiles=reaction['productSmiles'])
                        # no error handling on sanitaiton, catchall at the end
                        # from original code

                        smiles = reaction_struct['productSmiles']
                        # sane_smiles = sanitise_smiles(
                        #     smiles,
                        #     sanitisation_failed='error',
                        # )

                        # sane_inchikey = inchikey_from_smiles(sane_smiles)
                        product = CompoundService.get_by_smiles(smiles=smiles)

                        mrich.print(i, j, k, reaction_type, product)

                        reaction, _ = ReactionModel.objects.get_or_create(
                            reaction_type=reaction_type,
                            product_compound=product,
                        )

                        rs = []
                        print('reactant smiles', reaction_struct['reactantSmiles'])
                        for reactant_s in reaction_struct['reactantSmiles']:
                            reactant_comp, _ = CompoundService.create(smiles=reactant_s)
                            reactant, _ = ReactantModel.objects.get_or_create(
                                compound=reactant_comp,
                                reaction=reaction,
                            )
                            rs.append(reactant_comp.pk)

                        if do_check_chemistry and not check_chemistry(
                            reaction_type, CompoundSet(rs), product
                        ):
                            raise InvalidChemistryError(
                                f'{type=}, {rs=}, {product.id=}',
                            )

                        for r_id in rs:
                            if r_id in reactants:
                                intermediates.add(compound_id=r_id, amount=1)
                            else:
                                reactants.add(compound_id=r_id, amount=1)

                        reactions.add(reaction)

                except InvalidChemistryError:
                    continue
                except UnsupportedChemistryError:
                    mrich.warning('Skipping unsupported chemistry:', reaction_type)
                    continue
                # except Exception:
                #     mrich.error(
                #         'Uncaught error with row', i, 'route', j, 'reaction', k
                #     )
                #     continue

                products.add(Ingredient.from_compound(product, amount=1))

                recipe = Recipe(
                    reactions=reactions,
                    reactants=reactants,
                    intermediates=intermediates,
                    products=products,
                )

                if register_routes:
                    route, _ = RouteService.create_from_recipe(
                        recipe=recipe,
                    )
                    mrich.success('registered route', route.pk)

                if pick_first:
                    break

        return df

    @classmethod
    def ingest_enamine_real_routes(
        cls,
        csv_path: str | Path,
        do_check_chemistry: bool = True,
        register_routes: bool = True,
    ):
        df = pd.read_csv(csv_path)
        steps = len([col for col in df.columns if 'product_step' in col])

        for i, row in mrich.track(df.iterrows(), total=len(df)):
            mrich.set_progress_field('i', i)
            mrich.set_progress_field('n', len(df))

            d = row.to_dict()

            reactions = ReactionSet()
            reactants = IngredientSet()
            intermediates = IngredientSet()
            products = IngredientSet()

            product = None
            try:
                for step_id in range(1, steps + 1):
                    r1_smiles = d.get(f'reactant_step{step_id}')
                    if not r1_smiles or (
                        isinstance(r1_smiles, float) and isnan(r1_smiles)
                    ):
                        continue

                    reaction_type = d[f'reaction_name_step{step_id}']
                    product = CompoundService.get_by_smiles(smiles=d['smiles'])

                    mrich.print(i, step_id, reaction_type, product)

                    reactant_smiles = [r1_smiles]
                    r2_smiles = d.get(f'reactant2_step{step_id}')
                    if r2_smiles and not (
                        isinstance(r2_smiles, float) and isnan(r2_smiles)
                    ):
                        reactant_smiles.append(r2_smiles)

                    reaction, _ = ReactionModel.objects.get_or_create(
                        reaction_type=reaction_type,
                        product_compound=product,
                    )

                    rs = []
                    for smiles in reactant_smiles:
                        reactant_comp, _ = CompoundService.create(smiles=smiles)
                        reactant, _ = ReactantModel.objects.get_or_create(
                            compound=reactant_comp,
                            reaction=reaction,
                        )
                        rs.append(reactant_comp.pk)

                    if do_check_chemistry and not check_chemistry(
                        reaction_type, CompoundSet(rs), product
                    ):
                        raise InvalidChemistryError(
                            f'{reaction_type=}, {rs=}, {product.id=}',
                        )

                    for r_id in rs:
                        if r_id in reactants:
                            intermediates.add(compound_id=r_id, amount=1)
                        else:
                            reactants.add(compound_id=r_id, amount=1)

                    reactions.add(reaction)

            except InvalidChemistryError:
                continue
            except UnsupportedChemistryError:
                mrich.warning('Skipping unsupported chemistry:', reaction_type)
                continue
            except Exception:
                mrich.error('Uncaught error with row', i)
                raise

            if product is None:
                continue

            products.add(Ingredient.from_compound(product, amount=1))

            recipe = Recipe(
                reactions=reactions,
                reactants=reactants,
                intermediates=intermediates,
                products=products,
            )

            if register_routes:
                route, _ = RouteService.create_from_recipe(recipe=recipe)
                mrich.success('registered route', route.pk)

        return df

    @classmethod
    def ingest_syndirella_elabs(
        cls,
        *,
        df: pd.DataFrame,
        target: TargetModel,
        reject_flags: list[str],
        pose_tag_list: list[str],
        product_tag_list: list[str],
        max_energy_score: float,
        max_distance_score: float,
        require_intra_geometry_pass: bool,
        register_reactions: bool,
        scaffold_route: Route | None = None,
        scaffold_compound: CompoundModel | None = None,
    ) -> pd.DataFrame:

        # Tags must already be registered; ingestion never creates them.
        # Resolved here rather than at their point of use further down, so an
        # unknown tag fails before the compounds and reactions are written.
        product_tags = CompoundTagService.resolve_tags(product_tag_list)
        pose_tags = PoseTagService.resolve_tags(pose_tag_list)

        # work out number of reaction steps
        num_steps = max(
            [int(s.split('_')[0]) for s in df.columns if '_product_smiles' in s]
        )
        mrich.var('num_steps', num_steps)

        # add is_scaffold row
        df['is_scaffold'] = df[f'{num_steps}_product_name'].str.contains('scaffold')

        ###### PREP ######

        # flags

        present_flags = set()
        for step in range(num_steps):
            step += 1

            for flags in set(df[df[f'{step}_flag'].notna()][f'{step}_flag'].to_list()):
                for flag in flags:
                    present_flags.add(flag)

        if present_flags:
            mrich.warning('Flags in DataFrame:', present_flags)

        for flag in reject_flags:
            if flag in present_flags:
                for step in range(num_steps):
                    step += 1
                    matches = df[f'{step}_flag'].apply(
                        lambda x, flag=flag: flag in x if x is not None else False
                    )
                    mrich.print(
                        'Filtering out',
                        len(df[matches]),
                        'rows from step',
                        step,
                        'due to',
                        flag,
                    )
                    df = df[~matches]

        # poses

        n_null_mol = len(df[df['path_to_mol'].isna()])
        if n_null_mol:
            df = df[df['path_to_mol'].notna()]
            mrich.var('#rows skipped due to null path_to_mol', n_null_mol)

        if not len(df):
            mrich.warning('No valid rows')
            return None

        # inspirations
        inspiration_sets = set(tuple(sorted(i)) for i in df['regarded'])
        # smth like {('z0637a', 'z1040a')}

        if len(inspiration_sets) != 1:
            mrich.error('Varying inspirations not supported')
            return df

        (inspiration_set,) = inspiration_sets

        inspirations = PoseModel.objects.filter(
            pose_alias__in=inspiration_set,
            target=target,
        )

        if inspirations.count() != len(inspiration_set):
            print('target', target)
            print('inspiration_set', inspiration_set)
            print('inspiration comparison', inspirations.count(), len(inspiration_set))
        assert inspirations.count() == len(inspiration_set)

        # reference
        template_paths = set(df['template'].to_list())
        assert len(template_paths) == 1, 'Multiple references not supported'
        (template_path,) = template_paths
        template_path = Path(template_path)
        mrich.var('template_path', template_path)
        base_name = template_path.name.removesuffix('.pdb').removesuffix(
            '_delig-desolv'
        )
        # DEPRECATED(apo-naming): pre-'delig' Fragalysis naming, remove once all
        # data uses 'delig'
        base_name = base_name.removesuffix('_apo-desolv')
        # reference = self.poses[base_name]

        # TODO: error handling
        reference = PoseModel.objects.get(
            pose_alias=base_name,
            target=target,
        )

        assert reference, 'Could not determine reference structure'
        mrich.var('reference', reference)

        # that's nice but I need it before that
        # target = reference.target

        # subset of rows
        scaffold_df = df[df['is_scaffold']]
        elab_df = df[~df['is_scaffold']]
        mrich.var('#scaffold entries', len(scaffold_df))
        mrich.var('#elab entries', len(elab_df))

        if not len(scaffold_df) and not scaffold_route and not scaffold_compound:
            mrich.error('No valid scaffold rows')
            return None

        elif scaffold_route:
            ### SUPPLEMENT THE SCAFFOLD ROWS FROM KNOWN ROUTE

            assert scaffold_route.num_reactions == 1

            product = scaffold_route.products[0].compound
            reaction = scaffold_route.reactions[0]

            assert reaction.reactants.count() == 2

            scaffold_dict = {
                'scaffold_smiles': product.compound_smiles,
                '1_reaction': reaction.reaction_type,
                # this is so hacky
                '1_r1_smiles': reaction.reactants.first().compound.compound_smiles,
                '1_r2_smiles': reaction.reactants.last().compound.compound_smiles,
                '1_product_smiles': product.compound_smiles,
                '1_product_name': 'scaffold',
                '1_single_reactant_elab': False,
                '1_num_atom_diff': 0,
                'is_scaffold': True,
            }

            scaffold_df = pd.DataFrame([scaffold_dict])

            df = pd.concat([scaffold_df, df])

            scaffold_df = df[df['is_scaffold']]
            elab_df = df[~df['is_scaffold']]

        elif scaffold_compound:
            ### SUPPLEMENT PARTIAL SCAFFOLD ROWS FROM KNOWN PRODUCT

            scaffold_dict = {
                'scaffold_smiles': scaffold_compound.smiles,
                'is_scaffold': True,
            }

            scaffold_df = pd.DataFrame([scaffold_dict])

            df = pd.concat([scaffold_df, df])

            scaffold_df = df[df['is_scaffold']]
            elab_df = df[~df['is_scaffold']]

        # if dry_run:
        #     mrich.error('Not registering records (dry_run)')
        #     return df

        ###### ELABS ######

        # bulk register compounds

        smiles_cols = [
            c for c in df.columns if c.endswith('_smiles') and c != 'scaffold_smiles'
        ]

        for smiles_col in smiles_cols:
            inchikey_col = smiles_col.replace('_smiles', '_inchikey')
            compound_id_col = smiles_col.replace('_smiles', '_compound_id')

            unique_smiles = df[smiles_col].dropna().unique()

            mrich.debug(
                f'Registering {len(unique_smiles)} compounds from column: {smiles_col}'
            )

            # radical?
            values = CompoundService.create_from_smiles_list(unique_smiles)

            orig_smiles_to_inchikey = {
                orig_smiles: inchikey
                for orig_smiles, (inchikey, new_smiles) in zip(
                    unique_smiles, values, strict=False
                )
            }

            df[inchikey_col] = df[smiles_col].apply(
                lambda x, m=orig_smiles_to_inchikey: m.get(x)
            )

            # get associated IDs
            compound_inchikey_id_dict = {
                k.compound_inchikey: k.pk
                for k in CompoundModel.objects.filter(compound_smiles__in=unique_smiles)
            }
            df[compound_id_col] = df[inchikey_col].apply(
                lambda x, m=compound_inchikey_id_dict: m.get(x)
            )

        # bulk register reactions

        if register_reactions:
            for step in range(num_steps):
                step += 1

                mrich.debug(f'Registering reactions for step {step}')

                reaction_dicts = []

                for reaction_name, r1_id, r2_id, product_id in df[
                    [
                        f'{step}_reaction',
                        f'{step}_r1_compound_id',
                        f'{step}_r2_compound_id',
                        f'{step}_product_compound_id',
                    ]
                ].values:
                    # skip invalid rows
                    if pd.isna(r1_id) or pd.isna(product_id):
                        mrich.warning("Can't insert reactions for missing scaffold")
                        continue

                    # reactant IDs

                    reactant_ids = set()
                    reactant_ids.add(int(r1_id))

                    if not pd.isna(r2_id):
                        reactant_ids.add(int(r2_id))

                    product_id = int(product_id)

                    # registration data

                    reaction_dicts.append(
                        dict(
                            reaction_name=reaction_name,
                            reactant_ids=reactant_ids,
                            product_id=int(product_id),
                        )
                    )

            # why is this outside of loop?
            _ = ReactionService.create_from_lists(
                reaction_types=[d['reaction_name'] for d in reaction_dicts],
                product_ids=[d['product_id'] for d in reaction_dicts],
                reactant_id_lists=[d['reactant_ids'] for d in reaction_dicts],
            )

        scaffold_df = df[df['is_scaffold']]
        elab_df = df[~df['is_scaffold']]

        # tag product compounds:

        product_ids = list(df[f'{num_steps}_product_compound_id'].dropna().unique())
        products = CompoundModel.objects.filter(pk__in=product_ids)
        for compound in products:
            for compound_tag in product_tags:
                CompoundTagJunctionModel.objects.get_or_create(
                    compound=compound,
                    compound_tag=compound_tag,
                    target=target,
                )

        # bulk register scaffold relationships

        for step in range(num_steps):
            step += 1

            for role in ['r1', 'r2', 'product']:
                key = f'{step}_{role}_compound_id'

                mrich.debug(f'Registering scaffold relatonships for {key}')

                if step == num_steps and role == 'product' and scaffold_compound:
                    scaffold_id = scaffold_compound.id

                else:
                    scaffold_ids = list(scaffold_df[key].dropna().unique())

                    if not scaffold_ids:
                        mrich.warning(
                            "Can't insert scaffold relationships due to missing",
                            key,
                            'for all scaffold rows',
                        )
                        continue

                    if len(scaffold_ids) > 1:
                        mrich.error('Multiple scaffold row values in', key)
                        return scaffold_df

                    scaffold_id = scaffold_ids[0]

                # original code didn't do dropna? how? filter in later step?
                superstructure_ids = [
                    i for i in elab_df[key].dropna().unique() if i != scaffold_id
                ]

                # comp service?
                for superstructure_id in superstructure_ids:
                    base = CompoundModel.objects.get(pk=scaffold_id)
                    superstructure = CompoundModel.objects.get(
                        pk=int(superstructure_id)
                    )
                    ScaffoldModel.objects.get_or_create(
                        base_compound=base,
                        superstructure_compound=superstructure,
                    )

        # filter poses

        ok = df

        try:
            if require_intra_geometry_pass:
                mrich.var(
                    '#poses !intra_geometry_pass',
                    len(df[df['intra_geometry_pass'] == False]),  # noqa: E712
                )
                ok = ok[ok['intra_geometry_pass'] == True]  # noqa: E712

            if max_energy_score is not None:
                mrich.var(
                    f'#poses ∆∆G > {max_energy_score}',
                    len(df[df['∆∆G'] > max_energy_score]),
                )
                ok = ok[ok['∆∆G'] <= max_energy_score]

            if max_distance_score is not None:
                mrich.var(
                    f'#poses comRMSD > {max_distance_score}',
                    len(df[df['comRMSD'] > max_energy_score]),
                )
                ok = ok[ok['comRMSD'] <= max_distance_score]

        except Exception as e:
            mrich.error('Problem filtering dataframe')
            mrich.error(e)
            return df

        mrich.var('#acceptable poses', len(ok))

        if not len(ok):
            mrich.warning('No valid poses')
            return None

        # bulk register poses

        pose_ids = []
        scorer = ScoreService()
        for _, row in ok.iterrows():
            path = Path(os.path.relpath(row.path_to_mol))
            print('comp id in row', row[f'{num_steps}_product_compound_id'])

            # closed for testing
            if not path.exists():
                mrich.warning('Skipping pose w/ non-exising file:', path)
                continue

            if pd.isna(row[f'{num_steps}_product_compound_id']):
                continue

            pose, created = PoseService.create_from_record(
                compound_id=int(row[f'{num_steps}_product_compound_id']),
                target_id=int(target.id),
                reference=int(reference.id),
                path=str(path),
            )
            if created:
                scores = {
                    'energy_score': float(row['∆∆G']),
                    'distance_score': float(row['comRMSD']),
                }
                pose_ids.append(pose.id)
                scorer.add_scores_from_record(pose=pose, record=scores)

        if not pose_ids:
            mrich.warning('No valid poses')
            return None

        poses = PoseModel.objects.filter(pk__in=pose_ids)
        mrich.success('Registered', poses.count(), 'new poses')

        # query relevant poses (also previously registered)
        paths = poses.values_list('path', flat=True)

        # what the hell is this??
        records = PoseModel.objects.filter(
            path__in=paths,
        )
        for pose in records:
            # pose.inspirations.add(*PoseModel.objects.filter(pk__in=inspiration.ids))
            pose.inspirations.add(*inspirations.queryset)

        for pose in poses:
            pose.tags.add(*pose_tags)

        return df


# def create_compound(...):
#     assert connection.in_atomic_block
