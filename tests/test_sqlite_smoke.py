"""Smoke tests for the SQLite (container-free) test tier and its conftest fixtures.

Validates that hippo boots in SQLite mode and that compound registration works
offline -- i.e. that ``CompoundService.create`` populates the fields that the
Postgres cartridge triggers would otherwise fill.
"""

import pytest

pytestmark = pytest.mark.sqlite


def test_animal_boots_in_sqlite(animal):
    """The animal is created and has the expected target."""
    assert animal.target.target_name == "test"


def test_compounds_property_is_a_compound_set(animal):
    """animal.compounds returns a CompoundSet (empty or not)."""
    from designdb.sets.compound import CompoundSet

    assert isinstance(animal.compounds, CompoundSet)


def test_create_compound_populates_cartridge_fields(make_compound):
    """CompoundService.create fills mol/inchikey/hash in Python (no PG trigger)."""
    compound = make_compound("c1ccccc1O")  # phenol

    assert compound.pk is not None
    assert compound.compound_smiles == "c1ccccc1O"
    # populated in Python for the SQLite path
    assert compound.compound_inchikey == "ISWSIDIOOBJBQZ-UHFFFAOYSA-N"
    assert compound.compound_mol  # non-empty CTAB
    assert "RDKit" in compound.compound_mol
    assert compound.compound_hash  # tautomer-insensitive registration hash


def test_create_compound_is_idempotent(make_compound):
    """Re-registering the same SMILES returns the same compound (dedup by hash)."""
    a = make_compound("CCO")
    b = make_compound("CCO")

    assert a.pk == b.pk
