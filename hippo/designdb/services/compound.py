import logging
import re

import mrich
import rdkit
from designdb.models import CompoundModel, CompoundTagModel
from designdb.utils import (
    registration_hash_tautomer_insensitive,
    sanitise_smiles,
    superparent,
)
from django.conf import settings

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
    def tags_from_list(tag_list: list[str]):
        assert tag_list is not None, '"None" passed as tag_list'

        CompoundTagModel.objects.bulk_create(
            [
                CompoundTagModel(compound_tag_name=k.strip())
                for k in tag_list
                if k.strip()
            ],
            ignore_conflicts=True,
        )
        tags = CompoundTagModel.objects.filter(compound_tag_name__in=tag_list)
        return tags
