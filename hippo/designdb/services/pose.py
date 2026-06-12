import json
import logging
import re
from collections.abc import Iterable
from pathlib import Path

import mrich
import pandas as pd
import rdkit
# from rdkit.Chem import inchi
from designdb.models import CompoundModel, PoseModel, PoseTagModel, TargetModel
from designdb.utils import normalize_string_list
from designdb.utils_chem import get_rmsd
from designdb.utils_frag import GENERATED_TAG_COLS, META_IGNORE_COLS
from django.db.models import Q
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
        check_rmsd: bool = False,
        rmsd_threshold: float = 1.0,
    ):

        try:
            pose = PoseModel.objects.get(
                target=target,
                compound=compound,
                pose_alias=alias,
            )
            # default is to overwrite metadata. what about other props?
            # also, shoulnd't this be JSON?
            pose.metadata = metadata
            pose.save()
            created = False
        except PoseModel.DoesNotExist:
            if check_rmsd and (
                duplicate := cls.find_rmsd_duplicate(mol, compound, target, rmsd_threshold)
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
            created = True
        # except MultipleObjectsReturned:
        #     pass

        return pose, created

    @classmethod
    def find_rmsd_duplicate(
        cls,
        mol: 'Chem.rdchem.Mol',
        compound: 'CompoundModel',
        target: 'TargetModel',
        rmsd_threshold: float,
    ) -> 'PoseModel | None':
        for existing in PoseModel.objects.filter(compound=compound, target=target):
            try:
                rmsd = get_rmsd(mol, existing.pose_mol)
                if rmsd < rmsd_threshold:
                    logger.warning(
                        'Pose RMSD %.3f Å below threshold %.3f Å, skipping duplicate (alias=%s)',
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

    @staticmethod
    def tags_from_list(tag_list: list[str]):
        assert tag_list is not None, '"None" passed as tag_list'

        PoseTagModel.objects.bulk_create(
            [PoseTagModel(pose_tag_name=k.strip()) for k in tag_list if k.strip()],
            ignore_conflicts=True,
        )
        tags = PoseTagModel.objects.filter(pose_tag_name__in=tag_list)
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

        tags = PoseTagService.tags_from_list(pose_tag_set)

        return tags, metadata
