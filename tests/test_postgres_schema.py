"""DesignDB schema behaviour hippo relies on (Postgres tier).

These cover what the SQLite tier can't: the RDKit cartridge triggers that
populate compound/pose/catalogue fields, the hash-join triggers that link
catalogue prices to compounds (and so drive quoting), the score consistency
check, the score-pivot materialized view, and how a login without admin rights
is refused when registering tags/methods. Run with ``HIPPO_TEST_DB=postgres``
against a database from ``tests/postgres/start_db.sh``.
"""

from contextlib import contextmanager

import pytest

pytestmark = pytest.mark.postgres


def refetch(instance):
    """A fresh copy of ``instance`` from the DB (picks up trigger-set fields)."""
    return type(instance).objects.get(pk=instance.pk)


@pytest.fixture
def make_catalogue_price():
    """Factory: a catalogue entry with one price line, for a given SMILES/hash."""
    import rdkit
    from rdkit import Chem

    from designdb.models import CatalogueCompoundModel, CataloguePriceModel

    def _make(smiles: str, catalogue_hash: str, supplier: str = "Enamine"):
        catalogue = CatalogueCompoundModel.objects.create(
            catalogue_smiles=smiles,
            catalogue_hash=catalogue_hash,
            # NOT NULL in the schema, though nullable in the model
            rdkit_version=rdkit.__version__,
            inchi_version=Chem.inchi.GetInchiVersion(),
        )
        return CataloguePriceModel.objects.create(
            catalogue_compound=catalogue,
            vendor=supplier,
            supplier=supplier,
            supplier_id=f"{supplier}-{catalogue.pk}",
            amount=10,
            price=100,
            currency="USD",
        )

    return _make


### RDKit cartridge triggers


def test_compound_cartridge_fields_match_python(make_compound):
    """The cartridge fills compound_mol/inchikey on insert, and its InChIKey
    agrees with the one the SQLite tier computes in Python."""
    from rdkit import Chem
    from rdkit.Chem.inchi import MolToInchiKey

    smiles = "CC(=O)Nc1ccc(O)cc1"  # paracetamol
    compound = refetch(make_compound(smiles))

    assert compound.compound_inchikey == MolToInchiKey(Chem.MolFromSmiles(smiles))
    assert compound.compound_mol is not None


def test_pose_cartridge_fills_smiles_and_inchikey(animal, make_compound):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.inchi import MolToInchiKey

    from designdb.models import PoseModel

    smiles = "c1ccc2[nH]ccc2c1"  # indole
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    mol = Chem.RemoveHs(mol)

    pose = PoseModel.objects.create(
        compound=make_compound(smiles),
        target=animal.target,
        pose_alias="pg-pose-cartridge",
        pose_mol=mol,
    )
    pose = refetch(pose)

    assert Chem.CanonSmiles(pose.pose_smiles) == Chem.CanonSmiles(smiles)
    assert pose.pose_inchikey == MolToInchiKey(Chem.MolFromSmiles(smiles))


def test_catalogue_cartridge_rejects_invalid_smiles():
    from django.db import DatabaseError, transaction

    from designdb.models import CatalogueCompoundModel

    with pytest.raises(DatabaseError), transaction.atomic():
        CatalogueCompoundModel.objects.create(
            catalogue_smiles="not a smiles", catalogue_hash="pg-invalid"
        )


### catalogue <-> compound linking (drives quoting)


def test_price_added_after_compound_is_linked(make_compound, make_catalogue_price):
    from designdb.models import CompoundModel
    from designdb.sets.compound import CompoundSet

    compound = make_compound("CC(C)Cc1ccc(C(C)C(=O)O)cc1")  # ibuprofen
    make_catalogue_price("CC(C)Cc1ccc(C(C)C(=O)O)cc1", compound.compound_hash)

    cset = CompoundSet(CompoundModel.objects.filter(pk=compound.pk))

    assert list(cset.get_quoted().ids) == [compound.pk]
    assert len(cset.get_unquoted()) == 0


def test_compound_added_after_price_is_linked(make_compound, make_catalogue_price):
    from rdkit import Chem

    from designdb.services.quote import QuoteService
    from designdb.utils import registration_hash_tautomer_insensitive, superparent

    smiles = "OC(=O)c1ccccc1OC(C)=O"  # aspirin
    # the same hash CompoundService.create will compute, without registering
    h = registration_hash_tautomer_insensitive(superparent(Chem.MolFromSmiles(smiles)))

    price = make_catalogue_price(smiles, h)
    compound = make_compound(smiles)

    assert compound.compound_hash == h

    assert QuoteService.quoted_compound_ids([compound.pk]) == {compound.pk}
    assert compound.prices.filter(pk=price.pk).exists()


def test_unmatched_price_is_not_linked(make_compound, make_catalogue_price):
    from designdb.services.quote import QuoteService

    compound = make_compound("c1ccc2ccccc2c1")  # naphthalene
    make_catalogue_price("c1ccc2ccccc2c1", "pg-some-other-hash")

    quoted, unquoted = QuoteService.partition_quoted([compound.pk])

    assert quoted == set()
    assert unquoted == {compound.pk}


### scores


@pytest.fixture
def scored_pose(animal, make_compound):
    """A pose and a scoring method, for score-value tests."""
    from uuid import uuid4

    from designdb.models import PoseModel, ScoringMethodModel

    compound = make_compound("Cc1ccncc1")  # 4-methylpyridine
    pose = PoseModel.objects.create(
        compound=compound, target=animal.target, pose_alias="pg-scored"
    )
    # (method_name, method_version) is unique
    method = ScoringMethodModel.objects.create(
        method_name=f"pg_score_{uuid4().hex[:8]}", method_version="1.0"
    )
    return pose, method


def test_score_compound_must_match_pose(scored_pose, make_compound):
    from django.db import DatabaseError, transaction

    from designdb.models import ScoreValueModel

    pose, method = scored_pose
    other = make_compound("Cc1ccccn1")  # 2-methylpyridine

    with pytest.raises(DatabaseError), transaction.atomic():
        ScoreValueModel.objects.create(
            pose=pose, compound=other, scoring_method=method, score={"score": 1.0}
        )


def test_score_pivot_materialized_view(scored_pose):
    """Adding a scoring method adds a column to the pivot MV, and inserting a
    score refreshes it."""
    from django.db import connection

    from designdb.models import ScoreValueModel

    pose, method = scored_pose
    ScoreValueModel.objects.create(
        pose=pose,
        compound=pose.compound,
        scoring_method=method,
        score={"score": -7.5},
    )

    # column naming: see designdb.create_scores_per_pose_pivoted_mv()
    column = f"{method.method_name}_1_0_m{method.pk}"
    with connection.cursor() as cursor:
        cursor.execute(
            f"SELECT {column} FROM designdb.scores_per_pose_pivoted_mv "
            "WHERE pose_id = %s",
            [pose.pk],
        )
        rows = cursor.fetchall()

    assert len(rows) == 1
    assert float(rows[0][0]) == -7.5



### vocabulary permissions

READ_ONLY_ROLE = "hippo_test_read_only"


@contextmanager
def as_read_only():
    """Run as a role that may read, but not add to, the designdb tables.

    Stands in for a non-admin infrastructure login: the role is switched with
    ``SET ROLE`` on the test's own connection and reset afterwards.
    """
    from django.db import connection

    with connection.cursor() as cursor:
        cursor.execute("SELECT 1 FROM pg_roles WHERE rolname = %s", [READ_ONLY_ROLE])
        if cursor.fetchone() is None:
            cursor.execute(f"CREATE ROLE {READ_ONLY_ROLE} NOLOGIN")
        cursor.execute(f"GRANT USAGE ON SCHEMA designdb TO {READ_ONLY_ROLE}")
        cursor.execute(
            f"GRANT SELECT ON ALL TABLES IN SCHEMA designdb TO {READ_ONLY_ROLE}"
        )
        cursor.execute(f"SET ROLE {READ_ONLY_ROLE}")
    try:
        yield
    finally:
        with connection.cursor() as cursor:
            cursor.execute("RESET ROLE")


@pytest.mark.parametrize(
    "register",
    [
        lambda animal: animal.register_compound_tag("perm-new-compound-tag"),
        lambda animal: animal.register_pose_tag("perm-new-pose-tag"),
        lambda animal: animal.register_enumeration_method("perm-new-enum", "1.0"),
        lambda animal: animal.register_pose_method("perm-new-pose-method", "1.0"),
        lambda animal: animal.register_scoring_method("perm-new-scoring", "1.0"),
    ],
    ids=["compound_tag", "pose_tag", "enumeration", "pose_method", "scoring"],
)
def test_registration_without_admin_rights_is_explained(animal, register):
    from designdb.utils import InsufficientPrivilegeError

    with as_read_only():
        with pytest.raises(InsufficientPrivilegeError, match="DesignDB admin rights"):
            register(animal)


def test_registering_an_existing_tag_needs_no_admin_rights(animal):
    """get_or_create only reads when the tag exists, so non-admins can use it."""
    tag, created = animal.register_compound_tag("perm-existing-tag")
    assert created

    with as_read_only():
        again, created = animal.register_compound_tag("perm-existing-tag")

    assert not created
    assert again == tag
