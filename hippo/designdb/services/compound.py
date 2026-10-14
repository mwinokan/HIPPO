import logging
import re

import mrich
import rdkit
from designdb.models import CompoundModel, CompoundTagModel
from designdb.utils import (
    MissingTagError,
    admin_write,
    compound_hashes_from_smiles,
    registration_hash_tautomer_insensitive,
    safe_batch_size,
    sanitise_smiles,
    superparent,
)
from django.conf import settings
from django.db.models import QuerySet

# from mypackage.services.compound import CompoundService
from rdkit import Chem
from rdkit.Chem.inchi import MolToInchiKey

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


class CompoundBatchResult:
    def __init__(self):
        self.created = []
        self.errors = []


class CompoundService:
    @classmethod
    def create(
        cls,
        *,
        # mol: Chem.rdchem.Mol,
        smiles: str,
        # inchikey: str,
    ) -> tuple[CompoundModel, bool]:

        # designdb expects smils as input, so this is the entrypoint
        # for insertion
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        try:
            sp = superparent(mol)
        except Exception as e:
            raise ValueError(f'SuperParent failed: {e}') from e

        h = registration_hash_tautomer_insensitive(sp)

        defaults = {
            'compound_smiles': smiles,
            'rdkit_version': rdkit.__version__,
            'inchi_version': Chem.inchi.GetInchiVersion(),
        }

        # In SQLite mode there is no cartridge, so populate compound_mol (CTAB)
        # and compound_inchikey in Python. In Postgres the BEFORE INSERT trigger
        # (populate_compound_cartridge_from_smiles) fills these from the cartridge
        # and stays authoritative, so we leave them unset here.
        if settings.MANAGE_MODELS:
            defaults['compound_mol'] = Chem.MolToMolBlock(mol)
            defaults['compound_inchikey'] = MolToInchiKey(mol)

        compound, created = CompoundModel.objects.get_or_create(
            compound_hash=h,
            defaults=defaults,
        )
        if not created and logger.level == logging.DEBUG:
            mrich.warning(f'Skipping compound {h}, duplicate of {compound.pk}')

        # there's a following block in the original code
        # I don't understand what it is trying to achieve
        # smiles and inchikey are both inserted, so compound existing
        # but not reachable by inchikey should not happen. maybe this
        # covers compounds loaded through different pathway?

        # compound_id = self.db.insert_compound(
        #     smiles=smiles,
        #     tags=tags,
        #     warn_duplicate=debug,
        #     commit=False,
        # )

        # if not compound_id:
        #     inchikey = inchikey_from_smiles(smiles)
        #     compound = self.compounds[inchikey]

        #     if not compound:
        #         mrich.error(
        #             'CompoundModel exists in database but could not be found '
        #             'by inchikey'
        #         )
        #         mrich.var('smiles', smiles)
        #         mrich.var('inchikey', inchikey)
        #         mrich.var('observation_shortname', name)
        #         raise Exception

        # else:
        #     count_compound_registered += 1
        #     compound = self.compounds[compound_id]

        return compound, created

    @classmethod
    def create_batch(
        cls,
        *,
        smiles_list: list[str],
        max_workers: int | None = None,
        batch_size: int | None = None,
    ) -> dict[str, tuple[CompoundModel, bool]]:
        """Bulk counterpart of :meth:`create`.

        Resolves many SMILES to :class:`CompoundModel` rows in a fixed number of
        queries instead of one ``get_or_create`` (2 statements) per SMILES.

        The per-record version issues 2 statements per molecule; this issues at most
        4 in total regardless of batch size:

        0. ``SELECT`` existing compounds by exact SMILES
        1. ``SELECT`` existing compounds by registration hash
        2. ``bulk_create`` the missing ones (``ignore_conflicts``, so it is safe
           against a concurrent writer inserting the same hash)
        3. ``SELECT`` the full hash → compound mapping back

        Step 3 is needed because ``ignore_conflicts=True`` makes PostgreSQL skip the
        ``RETURNING`` clause, so ``bulk_create`` cannot populate primary keys. Steps
        1-3 are skipped entirely when step 0 resolves everything.

        Duplicate SMILES within ``smiles_list`` are collapsed before any DB work.

        Registration hashing dominates this method's cost (~6 ms/mol), so it is
        avoided wherever possible and parallelised otherwise:

        - ``compounds.compound_smiles`` is UNIQUE and indexed, so a SMILES already
          in the database is resolved by an indexed lookup and **never hashed**.
          Re-loading a file therefore costs almost no CPU.
        - whatever remains is hashed across processes by
          :func:`compound_hashes_from_smiles`.

        :param smiles_list: SMILES strings to resolve; duplicates are tolerated
        :param max_workers: hashing worker processes; see
            :func:`compound_hashes_from_smiles`
        :param batch_size: rows per INSERT statement; ``None`` picks the safe
            per-model maximum. See :func:`safe_batch_size`.
        :returns: mapping of SMILES to ``(compound, created)``, matching the return
            shape of :meth:`create`. SMILES whose hash could not be computed are
            absent from the mapping.
        """
        unique_smiles = list(dict.fromkeys(smiles_list))
        result: dict[str, tuple[CompoundModel, bool]] = {}

        # 0. Anything whose exact SMILES is already registered needs no hash at all.
        # compound_smiles carries a UNIQUE constraint, so this maps 1:1.
        known = CompoundModel.objects.filter(compound_smiles__in=unique_smiles)
        by_smiles = {c.compound_smiles: c for c in known}
        for smiles, compound in by_smiles.items():
            result[smiles] = (compound, False)

        unhashed = [s for s in unique_smiles if s not in by_smiles]
        if not unhashed:
            return result

        # SMILES -> hash for the remainder. A SMILES that cannot be hashed is
        # dropped rather than aborting the batch, mirroring how ingest_sdf skips
        # records it cannot sanitise.
        hash_by_smiles = compound_hashes_from_smiles(unhashed, max_workers=max_workers)
        for smiles in unhashed:
            if smiles not in hash_by_smiles:
                mrich.error(f'Could not hash {smiles=}')

        if not hash_by_smiles:
            return result

        wanted_hashes = set(hash_by_smiles.values())

        # 1. what already exists
        existing = {
            c.compound_hash: c
            for c in CompoundModel.objects.filter(compound_hash__in=wanted_hashes)
        }
        missing_hashes = wanted_hashes - set(existing)

        # 2. insert what does not
        if missing_hashes:
            # one representative SMILES per missing hash; distinct SMILES sharing a
            # hash are the same registered compound by definition
            smiles_for_hash = {}
            for smiles, h in hash_by_smiles.items():
                smiles_for_hash.setdefault(h, smiles)

            new_objects = []
            for h in missing_hashes:
                smiles = smiles_for_hash[h]
                defaults = {
                    'compound_hash': h,
                    'compound_smiles': smiles,
                    'rdkit_version': rdkit.__version__,
                    'inchi_version': Chem.inchi.GetInchiVersion(),
                }
                # In SQLite mode there is no cartridge, so populate compound_mol
                # (CTAB) and compound_inchikey in Python. In Postgres the BEFORE
                # INSERT trigger (populate_compound_cartridge_from_smiles) fills
                # these and stays authoritative, so we leave them unset here.
                if settings.MANAGE_MODELS:
                    mol = Chem.MolFromSmiles(smiles, sanitize=True)
                    defaults['compound_mol'] = Chem.MolToMolBlock(mol)
                    defaults['compound_inchikey'] = MolToInchiKey(mol)
                new_objects.append(CompoundModel(**defaults))

            CompoundModel.objects.bulk_create(
                new_objects,
                ignore_conflicts=True,
                batch_size=safe_batch_size(
                    CompoundModel, objs=new_objects, requested=batch_size
                ),
            )

        # 3. map every hash back to a row. Re-queried rather than reusing the
        # objects from step 2 because ignore_conflicts suppresses RETURNING.
        all_compounds = {
            c.compound_hash: c
            for c in CompoundModel.objects.filter(compound_hash__in=wanted_hashes)
        }

        for smiles, h in hash_by_smiles.items():
            compound = all_compounds.get(h)
            if compound is None:
                # only reachable if another transaction deleted the row between
                # steps 2 and 3
                mrich.error(f'Compound vanished during batch insert: {h}')
                continue
            result[smiles] = (compound, h in missing_hashes)

        return result

    # @classmethod
    # def create_from_smiles(
    #     cls,
    #     smiles: str,
    # ) -> tuple[CompoundModel, bool]:
    #     mol = Chem.MolFromSmiles(smiles, sanitize=True)
    #     compound, created = cls.create(mol=mol)
    #     return compound, created

    @classmethod
    def create_from_smiles_list(
        cls,
        smiles_list: list[str],
    ) -> list[tuple[str, str]]:
        result = []
        for smiles in smiles_list:
            sane_smiles = sanitise_smiles(
                smiles, verbosity=logger.level == logging.DEBUG
            )
            # compound, _ = cls.create_from_smiles(sane_smiles)
            compound, _ = cls.create(smiles=sane_smiles)
            result.append((compound.compound_inchikey, compound.compound_smiles))

        return result

    @classmethod
    def get_by_smiles(cls, smiles: str) -> CompoundModel | None:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        try:
            sp = superparent(mol)
        except Exception as e:
            raise ValueError(f'SuperParent failed: {e}') from e

        h = registration_hash_tautomer_insensitive(sp)

        return CompoundModel.objects.get(compound_hash=h)


class CompoundTagService:
    @staticmethod
    def register_tag(
        name: str, description: str | None = None, note: str | None = None
    ) -> tuple[CompoundTagModel, bool]:
        """Register a compound tag, or retrieve it if already registered.

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

        with admin_write(f'Registering compound tag {name!r}'):
            return CompoundTagModel.objects.get_or_create(
                compound_tag_name=name,
                defaults={
                    'compound_tag_description': description,
                    'compound_tag_note': note,
                },
            )

    @staticmethod
    def get_tags() -> QuerySet[CompoundTagModel]:
        """All registered compound tags."""
        return CompoundTagModel.objects.all()

    @staticmethod
    def resolve_tags(tag_list: list[str]) -> list[CompoundTagModel]:
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

        tags = list(CompoundTagModel.objects.filter(compound_tag_name__in=names))

        missing = names - {t.compound_tag_name for t in tags}
        if missing:
            raise MissingTagError(
                f'Unknown compound tag(s): {sorted(missing)}. Tags must exist '
                'before ingestion -- register them with '
                'animal.register_compound_tag() (needs DesignDB admin rights), '
                'then re-run.'
            )

        return tags
