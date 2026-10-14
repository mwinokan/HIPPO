import json
import logging
import re
from collections.abc import Iterable
from pathlib import Path

import mrich
import pandas as pd
import rdkit

# from rdkit.Chem import inchi
from designdb.models import (
    CompoundModel,
    PoseMethodJunctionModel,
    PoseMethodModel,
    PoseModel,
    PoseTagModel,
    TargetModel,
)
from designdb.utils import (
    MissingTagError,
    admin_write,
    normalize_string_list,
    safe_batch_size,
)
from designdb.utils_chem import get_rmsd
from designdb.utils_frag import GENERATED_TAG_COLS, META_IGNORE_COLS
from django.db.models import Q, QuerySet

# from mypackage.services.compound import CompoundService
from rdkit import Chem

# from .validation.compound import ValidationError, validate_compound_data

SDF_XCAv2_PATTERN = re.compile(
    r'^.*-.\d{4}_._\d*_\d_.*-.\d{4}\+.\+\d*\+\d_ligand\.sdf$'
)
SDF_XCAV3_PATTERN = re.compile(
    r'^.*-.\d{4}_._\d*_._\d_.*-.\d{4}\+.\+\d*\+.\+\d_ligand\.sdf$'
)


SDF_FRAGALYSIS_PATTERN = re.compile(r'^.*\d{4}[a-z].sdf$')
PDBID_PATTERN = re.compile(r'^[A-Za-z0-9]{4}-[a-z].sdf$')


logger = logging.getLogger(__name__)


class PoseService:
    @classmethod
    def create(
        cls,
        *,
        compound: CompoundModel,
        target: TargetModel,
        mol: Chem.rdchem.Mol,
        alias: str,
        path: str,
        metadata: dict[str, str],
        inchikey: str,
        smiles: str,
        reference: int | None = None,
        pose_method: 'PoseMethodModel | None' = None,
        check_rmsd: bool = False,
        rmsd_threshold: float = 1.0,
    ):
        qs = PoseModel.objects.filter(
            target=target,
            compound=compound,
            pose_alias=alias,
        )
        if pose_method:
            qs = qs.filter(methods=pose_method)

        # (target, compound, pose_alias, method) should identify one pose
        try:
            existing = qs.get()
        except PoseModel.DoesNotExist:
            existing = None
        else:
            # overwrite metadata on the matching pose
            existing.pose_metadata = json.dumps(metadata)
            existing.save()
            return existing, False

        if check_rmsd and (
            duplicate := cls.find_rmsd_duplicate(
                mol, compound, target, rmsd_threshold, pose_method=pose_method
            )
        ):
            return duplicate, False

        pose = PoseModel(
            compound=compound,
            target=target,
            pose_alias=alias,
            protein_link=path,
            pose_inchikey=inchikey,  # SQLITE_RELIC
            pose_smiles=smiles,  # SQLITE_RELIC
            pose_metadata=json.dumps(metadata),
            pose_mol=mol,
            # pose_mol=Chem.MolToMolBlock(mol),
            rdkit_version=rdkit.__version__,
            inchi_version=Chem.inchi.GetInchiVersion(),
            pose_reference=reference,
        )
        pose.save()

        # associate the method here so subsequent loads dedup correctly
        if pose_method is not None:
            pose.methods.add(pose_method)

        return pose, True

    @staticmethod
    def resolve_aliases_batch(
        aliases: Iterable[str],
        target: TargetModel | None = None,
    ) -> dict[str, PoseModel]:
        """Resolve many pose aliases to poses in a single query.

        Shared lookup behind :meth:`get_inspirations_batch` and
        :meth:`get_reference_batch`. The per-record versions issue one query per
        record; this issues one for the whole batch.

        .. note::
           Until an index exists on ``poses (pose_alias, target_id)`` this query is
           a sequential scan -- but it is now one scan per batch rather than one per
           record. See ``docs/proposal_pose_alias_index.md``.

        :param aliases: pose aliases to look up
        :param target: optional target to scope the lookup to
        :returns: mapping of alias to pose; aliases with no match are absent
        """
        aliases = {a for a in aliases if a}
        if not aliases:
            return {}

        qs = PoseModel.objects.filter(pose_alias__in=aliases)
        if target is not None:
            qs = qs.filter(target=target)

        # last row wins on duplicate aliases, matching the arbitrary pick that
        # .get() would make (it would actually raise; see get_reference_batch)
        return {p.pose_alias: p for p in qs}

    @classmethod
    def get_inspirations_batch(
        cls,
        records: list[dict],
        *,
        global_inspirations: list[int],
        inspiration_map: dict,
        inspiration_col: str | None,
        name_col: str,
        target: TargetModel | None = None,
    ) -> dict[int, list[int]]:
        """Bulk counterpart of :meth:`get_inspirations`.

        The per-record version runs one ``pk__in ... OR pose_alias__in ...`` query
        per record. That ``OR`` cannot use a bitmap index scan while ``pose_alias``
        is unindexed, so each record pays a full scan of ``poses``. This collects
        every pk and alias across the whole batch first, then resolves them in at
        most two queries total.

        :param records: preprocessed SDF records
        :param global_inspirations: inspiration ids applied to every record
        :param inspiration_map: mapping of record name to inspiration pose list
        :param inspiration_col: record column holding per-record inspirations
        :param name_col: record column holding the record name
        :param target: target to scope alias lookups to
        :returns: mapping of record index to a list of inspiration pose ids
        """
        # pass 1: parse every record's inspiration references, no DB access
        per_record_refs: dict[int, list] = {}
        all_pks: set[int] = set()
        all_aliases: set[str] = set()

        for idx, r in enumerate(records):
            parsed = []
            sources = [
                global_inspirations,
                inspiration_map.get(r[name_col], []),
            ]
            if inspiration_col:
                sources.append(r.get(inspiration_col, []))

            for el in sources:
                if el is None:
                    continue
                if isinstance(el, str):
                    parsed.extend(normalize_string_list(el))
                elif isinstance(el, Iterable) and not isinstance(el, dict):
                    parsed.extend(el)
                else:
                    logger.warning(
                        'Unsupported inspiration collection received: %s',
                        type(el),
                    )

            refs = []
            for val in parsed:
                try:
                    pk = int(val)
                except (TypeError, ValueError):
                    refs.append(('alias', val))
                    all_aliases.add(val)
                else:
                    refs.append(('pk', pk))
                    all_pks.add(pk)
            per_record_refs[idx] = refs

        # pass 2: resolve all references in bulk
        valid_pks = (
            set(PoseModel.objects.filter(pk__in=all_pks).values_list('pk', flat=True))
            if all_pks
            else set()
        )
        by_alias = cls.resolve_aliases_batch(all_aliases, target=target)

        # pass 3: reassemble per record, preserving order and dropping duplicates
        result: dict[int, list[int]] = {}
        for idx, refs in per_record_refs.items():
            ids: list[int] = []
            seen: set[int] = set()
            for kind, val in refs:
                if kind == 'pk':
                    pose_id = val if val in valid_pks else None
                else:
                    pose = by_alias.get(val)
                    pose_id = pose.pk if pose else None
                if pose_id is not None and pose_id not in seen:
                    seen.add(pose_id)
                    ids.append(pose_id)
            result[idx] = ids

        return result

    @classmethod
    def get_reference_batch(
        cls,
        records: list[dict],
        *,
        reference_col: str,
        target: TargetModel,
    ) -> dict[int, int | None]:
        """Bulk counterpart of :meth:`get_reference`.

        .. warning::
           This fixes a bug in the per-record path. In ``ingest_sdf`` the reference
           is resolved with ``if not reference and reference_col:``, assigning to the
           *function parameter*. Once the first record sets it, the condition is
           false forever after, so **every** pose in the file silently inherits the
           first record's reference. This version resolves each record's own
           reference, so batch results will legitimately differ from the per-record
           path for any SDF whose ``reference_col`` is not constant.

        :param records: preprocessed SDF records
        :param reference_col: record column holding the reference alias or id
        :param target: target to scope alias lookups to
        :returns: mapping of record index to reference pose id (``None`` if absent)
        """
        raw_by_idx: dict[int, str | int] = {}
        aliases: set[str] = set()
        pks: set[int] = set()

        for idx, r in enumerate(records):
            raw = r.get(reference_col)
            if raw is None or (isinstance(raw, float) and pd.isna(raw)):
                continue
            raw_by_idx[idx] = raw
            try:
                pks.add(int(raw))
            except (TypeError, ValueError):
                aliases.add(str(raw))

        by_alias = cls.resolve_aliases_batch(aliases, target=target)
        valid_pks = (
            set(PoseModel.objects.filter(pk__in=pks).values_list('pk', flat=True))
            if pks
            else set()
        )

        result: dict[int, int | None] = {}
        for idx, raw in raw_by_idx.items():
            try:
                pk = int(raw)
            except (TypeError, ValueError):
                pose = by_alias.get(str(raw))
                result[idx] = pose.pk if pose else None
            else:
                result[idx] = pk if pk in valid_pks else None

        return result

    @classmethod
    def create_batch(
        cls,
        *,
        target: TargetModel,
        specs: list[dict],
        pose_method: PoseMethodModel | None = None,
        batch_size: int | None = None,
        check_rmsd: bool = False,
        rmsd_threshold: float = 1.0,
    ) -> tuple[list[PoseModel | None], int]:
        """Bulk counterpart of :meth:`create`.

        The per-record version issues a lookup, an insert and (for pose methods) a
        junction insert per pose. This resolves all existing poses in one query,
        inserts the new ones in one ``bulk_create``, and attaches methods in one
        more.

        Each spec is a dict with the same keys :meth:`create` takes as kwargs:
        ``compound``, ``mol``, ``alias``, ``path``, ``metadata``, ``inchikey``,
        ``smiles``, ``reference``.

        Matches :meth:`create` semantics: an existing pose has its metadata
        overwritten and is returned with ``created=False``. Intra-batch duplicates
        (two records with the same compound and alias) collapse to one pose, with
        the last record's metadata winning -- the same end state the sequential
        version reaches.

        ``check_rmsd`` is supported. The pairwise RMSD comparison is unavoidably
        per-pose CPU work, but the *candidate lookup* batches like everything else:
        one query fetches every existing pose for the batch's compounds, and each
        spec is then compared against its own compound's candidates in memory.
        Sequential semantics are preserved -- :meth:`create` saves each pose before
        testing the next, so a pose can match one created earlier in the same file;
        specs accepted earlier in this batch are therefore included as candidates
        for later ones.

        :param target: target the poses belong to
        :param specs: per-pose field dicts, in record order
        :param pose_method: optional method to associate with every created pose
        :param batch_size: rows per INSERT/UPDATE statement; ``None`` picks the safe
            per-model maximum. See :func:`safe_batch_size`.
        :param check_rmsd: skip a pose whose RMSD to an existing pose of the same
            compound is below ``rmsd_threshold``
        :param rmsd_threshold: RMSD below which two poses are the same, in Angstrom
        :returns: ``(poses, created_count)`` where ``poses`` is aligned with
            ``specs`` and holds ``None`` for any entry that could not be created
        """
        if not specs:
            return [], 0

        # collapse intra-batch duplicates; last spec for a key wins
        key_for_idx: list[tuple[int, str]] = [
            (spec['compound'].pk, spec['alias']) for spec in specs
        ]
        spec_for_key: dict[tuple[int, str], dict] = {}
        for key, spec in zip(key_for_idx, specs, strict=True):
            spec_for_key[key] = spec

        compound_ids = {k[0] for k in spec_for_key}
        aliases = {k[1] for k in spec_for_key}

        # 1. one query for everything that already exists
        existing_qs = PoseModel.objects.filter(
            target=target,
            compound_id__in=compound_ids,
            pose_alias__in=aliases,
        )
        if pose_method is not None:
            existing_qs = existing_qs.filter(methods=pose_method)

        existing: dict[tuple[int, str], PoseModel] = {}
        for pose in existing_qs:
            existing[(pose.compound_id, pose.pose_alias)] = pose

        # 2. overwrite metadata on the ones that already exist
        to_update = []
        for key, pose in existing.items():
            spec = spec_for_key.get(key)
            if spec is None:
                continue
            pose.pose_metadata = json.dumps(spec['metadata'])
            to_update.append(pose)
        if to_update:
            PoseModel.objects.bulk_update(
                to_update,
                ['pose_metadata'],
                batch_size=safe_batch_size(
                    PoseModel,
                    objs=to_update,
                    fields=['pose_metadata'],
                    requested=batch_size,
                ),
            )

        # 3. insert the ones that do not
        new_keys = [k for k in spec_for_key if k not in existing]

        # RMSD dedup, for specs that did not match on (compound, alias). One query
        # for all candidates, then per-spec comparison. `rmsd_match` maps a skipped
        # key either to an existing pose or to another key in this batch whose pose
        # is created below.
        rmsd_match: dict[tuple[int, str], PoseModel] = {}
        rmsd_alias: dict[tuple[int, str], tuple[int, str]] = {}
        if check_rmsd and new_keys:
            cand_qs = PoseModel.objects.filter(
                target=target, compound_id__in={k[0] for k in new_keys}
            )
            if pose_method is not None:
                cand_qs = cand_qs.filter(methods=pose_method)

            candidates: dict[int, list] = {}
            for pose in cand_qs:
                candidates.setdefault(pose.compound_id, []).append(
                    (pose, pose.pose_mol)
                )

            remaining = []
            for key in new_keys:
                spec = spec_for_key[key]
                mol = spec['mol']
                hit = None
                for pose, other in candidates.get(key[0], []):
                    try:
                        if get_rmsd(mol, other) < rmsd_threshold:
                            hit = pose
                            break
                    except Exception:
                        logger.warning('RMSD calculation failed against pose %s', pose)
                if hit is not None:
                    if isinstance(hit, PoseModel):
                        logger.warning(
                            'Pose RMSD below threshold %.3f A, skipping duplicate '
                            '(alias=%s)',
                            rmsd_threshold,
                            spec['alias'],
                        )
                        rmsd_match[key] = hit
                    else:
                        rmsd_alias[key] = hit
                    continue
                # accepted: becomes a candidate for later specs, mirroring the
                # sequential path where each pose is saved before the next is tested
                candidates.setdefault(key[0], []).append((key, mol))
                remaining.append(key)
            new_keys = remaining

        created_poses: dict[tuple[int, str], PoseModel] = {}
        if new_keys:
            new_objects = []
            for key in new_keys:
                spec = spec_for_key[key]
                new_objects.append(
                    PoseModel(
                        compound=spec['compound'],
                        target=target,
                        pose_alias=spec['alias'],
                        protein_link=spec['path'],
                        pose_inchikey=spec['inchikey'],  # SQLITE_RELIC
                        pose_smiles=spec['smiles'],  # SQLITE_RELIC
                        pose_metadata=json.dumps(spec['metadata']),
                        pose_mol=spec['mol'],
                        rdkit_version=rdkit.__version__,
                        inchi_version=Chem.inchi.GetInchiVersion(),
                        pose_reference=spec.get('reference'),
                    )
                )
            # no ignore_conflicts: PostgreSQL only populates primary keys via
            # RETURNING, which is suppressed when conflicts are ignored, and step 1
            # has already excluded the rows that exist
            inserted = PoseModel.objects.bulk_create(
                new_objects,
                batch_size=safe_batch_size(
                    PoseModel, objs=new_objects, requested=batch_size
                ),
            )
            for key, pose in zip(new_keys, inserted, strict=True):
                created_poses[key] = pose

            # 4. attach the method to every new pose in one statement
            if pose_method is not None:
                PoseMethodJunctionModel.objects.bulk_create(
                    [
                        PoseMethodJunctionModel(pose=pose, pose_method=pose_method)
                        for pose in created_poses.values()
                    ],
                    ignore_conflicts=True,
                    batch_size=safe_batch_size(
                        PoseMethodJunctionModel, requested=batch_size
                    ),
                )

        resolved = {**existing, **created_poses, **rmsd_match}
        # keys that matched a pose created within this same batch
        for key, other_key in rmsd_alias.items():
            if other_key in resolved:
                resolved[key] = resolved[other_key]
        poses = [resolved.get(key) for key in key_for_idx]
        return poses, len(created_poses)

    @classmethod
    def find_rmsd_duplicate(
        cls,
        mol: 'Chem.rdchem.Mol',
        compound: 'CompoundModel',
        target: 'TargetModel',
        rmsd_threshold: float,
        pose_method: 'PoseMethodModel | None' = None,
    ) -> 'PoseModel | None':
        # only compare against poses produced by the same method
        candidates = PoseModel.objects.filter(compound=compound, target=target)
        if pose_method is not None:
            candidates = candidates.filter(methods=pose_method)
        for existing in candidates:
            try:
                rmsd = get_rmsd(mol, existing.pose_mol)
                if rmsd < rmsd_threshold:
                    logger.warning(
                        'Pose RMSD %.3f Å below threshold %.3f Å, '
                        'skipping duplicate (alias=%s)',
                        rmsd,
                        rmsd_threshold,
                        existing.pose_alias,
                    )
                    return existing
            except Exception:
                logger.warning('RMSD calculation failed for pose pk=%s', existing.pk)
        return None

    @classmethod
    def create_from_record(
        cls,
        *,
        compound_id: int,
        target_id: int,
        path: str,
        reference: int | None = None,
    ):
        target = TargetModel.objects.get(pk=target_id)
        compound = CompoundModel.objects.get(pk=compound_id)
        pose, created = PoseModel.objects.get_or_create(
            compound=compound,
            target=target,
            protein_link=path,
            reference=reference,
        )
        return pose, created

    # this is parsing input, maybe in ingestion?
    @staticmethod
    def get_inspirations(*args, target: TargetModel | None = None):
        parsed = []
        for el in args:
            if isinstance(el, str):
                parsed.extend(normalize_string_list(el))
            elif isinstance(el, Iterable) and not isinstance(el, dict):
                parsed.extend(el)
            else:
                logger.warning(
                    'Unsupported inspiration collection received: %s',
                    type(el),
                )

        # inputs can be pk or name
        pks = []
        aliases = []

        for val in parsed:
            try:
                pks.append(int(val))
            except ValueError:
                # assume string alias
                aliases.append(val)

        qs = PoseModel.objects.filter(
            Q(pk__in=pks) | Q(pose_alias__in=aliases, target=target)
        )

        return qs

    @staticmethod
    def get_reference(reference, target) -> int:
        try:
            reference = int(reference)
            # should I check if exist here as well?
        except ValueError:
            try:
                reference = PoseModel.objects.get(
                    pose_alias=reference,
                    target=target,
                ).pk
            except PoseModel.DoesNotExist as exp:
                logger.error('PoseModel %s does not exist', reference)
                raise PoseModel.DoesNotExist from exp

        return reference


class PoseTagService:
    def __init__(self, metadata_file: Path | str, other_tags: list[str] | None = None):
        self._df = pd.read_csv(metadata_file)
        self._curated_tag_cols = [
            c
            for c in self._df.columns
            if c not in META_IGNORE_COLS + GENERATED_TAG_COLS
        ]
        # any other tags to be added
        if other_tags:
            self._other_tags = [k.strip() for k in other_tags if k.strip()]
        else:
            self._other_tags = []

        mrich.var('curated_tag_cols', self._curated_tag_cols)

        # Resolve the whole vocabulary once, up front: an unregistered tag fails
        # before any record is written rather than part-way through the load, and
        # tags_and_meta stops issuing a query per record. Only columns that are
        # truthy somewhere are resolved -- an all-false column never produced a
        # tag, so requiring it to exist would be stricter than the old behaviour.
        # NB the truthiness test must match the one tags_and_meta applies per
        # record below -- pandas' Series.any() skips NaN, but bool(nan) is True,
        # so an all-NaN column would be left out of the map and then asked for.
        used_tag_cols = [
            c
            for c in self._curated_tag_cols
            if any(bool(v) for v in self._df[c].values)
        ]
        self._tag_map = {
            tag.pose_tag_name: tag
            for tag in PoseTagService.resolve_tags(used_tag_cols + self._other_tags)
        }

    @staticmethod
    def register_tag(
        name: str, description: str | None = None, note: str | None = None
    ) -> tuple[PoseTagModel, bool]:
        """Register a pose tag, or retrieve it if already registered.

        Registration counterpart of :meth:`resolve_tags`, in the same spirit as
        :meth:`.MethodService.register_pose_method`. Adding to the vocabulary
        needs DesignDB admin rights; retrieving an existing tag does not.

        :param name: tag name; surrounding whitespace is ignored
        :param description: description, used only when the tag is created
        :param note: note, used only when the tag is created
        :returns: ``(tag, created)``
        :raises ValueError: if ``name`` is blank
        :raises InsufficientPrivilegeError: if the tag is new and the database
            login may not add it
        """
        name = name.strip() if name else ''
        if not name:
            raise ValueError('Tag name must not be blank')

        with admin_write(f'Registering pose tag {name!r}'):
            return PoseTagModel.objects.get_or_create(
                pose_tag_name=name,
                defaults={
                    'pose_tag_description': description,
                    'pose_tag_note': note,
                },
            )

    @staticmethod
    def get_tags() -> QuerySet[PoseTagModel]:
        """All registered pose tags."""
        return PoseTagModel.objects.all()

    @staticmethod
    def resolve_tags(tag_list: list[str]) -> list[PoseTagModel]:
        """Resolve tag names to existing rows in the tag vocabulary.

        Lookup-only: an unknown name is an error rather than a new row; tags are
        added with :meth:`register_tag`, in the same spirit as
        :meth:`.MethodService.resolve_pose_method`.

        :param tag_list: tag names; surrounding whitespace, blanks and duplicates
            are ignored
        :returns: the matching tags, evaluated (not a lazy queryset -- callers
            iterate them once per chunk)
        :raises MissingTagError: if any name is not in the vocabulary
        """
        if tag_list is None:
            raise ValueError('"None" passed as tag_list')

        names = {k.strip() for k in tag_list if k and k.strip()}
        if not names:
            return []

        tags = list(PoseTagModel.objects.filter(pose_tag_name__in=names))

        missing = names - {t.pose_tag_name for t in tags}
        if missing:
            raise MissingTagError(
                f'Unknown pose tag(s): {sorted(missing)}. Tags must exist before '
                'ingestion -- register them with animal.register_pose_tag() '
                '(needs DesignDB admin rights), then re-run.'
            )

        return tags

    # might be a good idea to break meta and tags apart
    def tags_and_meta(
        self,
        *,
        code: str,
        longcode: str,
    ) -> tuple[list[PoseTagModel], dict[str, str]]:
        meta_row = self._df[self._df['Code'] == code]
        if not len(meta_row):
            meta_row = self._df[self._df['Long code'] == longcode]

        # TODO: another unhandled exception, apprently not having
        # meta_row is an option

        metadata = {'fragalysis_longcode': meta_row['Long code'].values[0]}

        for tag in GENERATED_TAG_COLS:
            if tag in meta_row.columns:
                metadata[tag] = meta_row[tag].values[0]

        pose_tag_set = set(self._other_tags)

        for tag in self._curated_tag_cols:
            if meta_row[tag].values[0]:
                pose_tag_set.add(tag)

        # keys come from resolve_tags, which strips; column names may not be
        tags = [self._tag_map[tag.strip()] for tag in pose_tag_set]

        return tags, metadata
