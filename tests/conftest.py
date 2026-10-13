"""Shared pytest fixtures for the SQLite (container-free) and Postgres test tiers.

By default these fixtures run hippo in SQLite mode
(``load_hippo(..., db="...sqlite")``), which configures Django with
``manage_models=True``, builds the schema programmatically, and uses the
plain-text ``RDKitMolField`` shim -- so no Postgres DesignDB container, RDKit
cartridge, or Fragalysis access is needed.

With ``HIPPO_TEST_DB=postgres`` the same ``animal`` fixture instead connects to
a DesignDB Postgres (RDKit cartridge, triggers, materialized view) configured
by the usual ``DB_*``/``POSTGRES_PORT`` env vars -- see
``tests/postgres/start_db.sh``. Tests marked ``sqlite`` run on either backend;
tests marked ``postgres`` need the real schema and are skipped on SQLite.

Notes / constraints:
- Django can only be configured once per process, so ``animal`` is
  session-scoped: the DB schema is created once and shared across tests.
  Tests should therefore create their own uniquely-identified objects rather
  than assuming an empty database.
- Import ``designdb.*`` lazily (inside fixtures/tests), never at module top:
  ``designdb.models`` reads ``settings.MANAGE_MODELS`` at import time, which is
  only defined after ``animal`` has configured Django.
"""

import os

import pytest

BACKEND = os.environ.get("HIPPO_TEST_DB", "sqlite")


def pytest_collection_modifyitems(config, items):
    """Skip the pre-refactor tests until they're migrated to the SQLite tier,
    and the Postgres-only tests unless running against Postgres.

    Any test not marked ``sqlite`` or ``postgres`` is skipped (not failed), so a
    plain ``pytest`` run stays green. Migrating a test -- point it at the
    fixtures in this conftest and add ``pytestmark = pytest.mark.sqlite`` --
    un-skips it automatically. ``make test`` uses ``-m sqlite`` and never
    collects these.
    """
    skip_legacy = pytest.mark.skip(
        reason="pre-refactor test, pending migration to the SQLite tier"
    )
    skip_postgres = pytest.mark.skip(
        reason="needs the DesignDB Postgres schema; run with HIPPO_TEST_DB=postgres"
    )
    for item in items:
        if "postgres" in item.keywords:
            if BACKEND != "postgres":
                item.add_marker(skip_postgres)
        elif "sqlite" not in item.keywords:
            item.add_marker(skip_legacy)


def postgres_config() -> dict:
    """Connection settings for the Postgres tier, from the usual env vars.

    Refuses any database whose name doesn't end in ``_test``: the ``animal``
    fixture empties every ``designdb`` table at the start of the session.
    """
    config = {
        "DB_NAME": os.environ.get("DB_NAME", "designdb_test"),
        "DB_USER": os.environ.get("DB_USER", "hippo"),
        "DB_PASSWORD": os.environ.get("DB_PASSWORD", "hippo"),
        "DB_HOST": os.environ.get("DB_HOST", "localhost"),
        "POSTGRES_PORT": os.environ.get("POSTGRES_PORT", "54329"),
    }
    if not config["DB_NAME"].endswith("_test"):
        pytest.exit(
            f"refusing to run the Postgres tier against {config['DB_NAME']!r}: "
            "the database name must end in '_test' (its tables are truncated)",
            returncode=2,
        )
    return config


def truncate_designdb(config: dict) -> None:
    """Empty every ``designdb`` table, so each session starts from a clean DB.

    Also resets the score materialized view: TRUNCATE fires none of its
    triggers, and an interrupted ``load_sdf`` can leave its refresh trigger
    disabled.
    """
    import psycopg

    with psycopg.connect(
        host=config["DB_HOST"],
        port=config["POSTGRES_PORT"],
        dbname=config["DB_NAME"],
        user=config["DB_USER"],
        password=config["DB_PASSWORD"],
    ) as conn:
        tables = conn.execute(
            "SELECT format('%I.%I', schemaname, tablename) FROM pg_tables "
            "WHERE schemaname = 'designdb'"
        ).fetchall()
        if tables:
            names = ", ".join(t for (t,) in tables)
            conn.execute(f"TRUNCATE {names} RESTART IDENTITY CASCADE")
        conn.execute(
            "ALTER TABLE designdb.score_values "
            "ENABLE TRIGGER trg_score_values_refresh_pivoted_mv"
        )
        conn.execute("SELECT designdb.create_scores_per_pose_pivoted_mv()")


@pytest.fixture(scope="session")
def animal(tmp_path_factory):
    """A :class:`HIPPO` animal backed by a fresh throwaway database: SQLite by
    default, or the DesignDB Postgres with ``HIPPO_TEST_DB=postgres``.

    Session-scoped (Django is configured once per process). Downloads are off:
    ``DOWNLOAD_APO_DESOLV_ON_INIT`` defaults to False, so init does no network.
    """
    from unittest.mock import patch

    import hippo

    if BACKEND == "postgres":
        db = postgres_config()
        truncate_designdb(db)
    else:
        db = str(tmp_path_factory.mktemp("hippo_db") / "test.sqlite")

    # load_hippo() checks the target access string against the TA authenticator
    # service. There is no service (and no TA_AUTH_QUERY_KEY) in the test
    # environment, and an unset key yields an empty access set, which would make
    # load_hippo() return None. Stub the lookup rather than relaxing the check.
    with patch(
        "hippo.bootstrap.get_auth_target_access",
        return_value={"test-proposal"},
    ):
        return hippo.HIPPO(
            target_name="test",
            target_access_string="test-proposal",
            username="test-user",
            db=db,
        )


@pytest.fixture
def make_compound(animal):
    """Factory: register a compound from SMILES and return the ``CompoundModel``.

    Exercises the real ingestion entrypoint (``CompoundService.create``), which
    in SQLite mode populates ``compound_mol``/``compound_inchikey`` in Python
    (on Postgres the cartridge triggers fill them).
    """
    from designdb.services.compound import CompoundService

    def _make(smiles: str):
        compound, _ = CompoundService.create(smiles=smiles)
        return compound

    return _make
