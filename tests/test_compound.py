"""Compound component properties, modernized onto the SQLite ``animal`` fixture.

Replaces the pre-refactor version (positional ``hippo.HIPPO('test', DB)``,
``animal.C1``, legacy ``.db`` property, ``animal.db.close()``). The compound is
registered through the real ingestion entrypoint (``CompoundService.create``,
via the ``make_compound`` fixture) and wrapped in the new ``Compound`` component.

The fixture compound has no poses/reactions/scaffolds, so properties that need
related data (e.g. ``best_placed_pose``) are covered by data-backed tests
elsewhere, not here.
"""

import pytest

pytestmark = pytest.mark.sqlite

# Must be populated for a freshly-registered compound (no related data).
# 0 / False / empty-collection all count as "not None".
NOT_NULL_PROPERTIES = [
    "id",
    "inchikey",
    "name",
    "smiles",
    "mol",
    "num_heavy_atoms",
    "molecular_weight",
    "num_rings",
    "formula",
    "atomtype_dict",
    "tags",
    "poses",
    "reactions",
    "num_poses",
    "num_reactions",
    "num_reactant",
    "num_scaffolds",
    "is_scaffold",
    "is_elab",
    "is_product",
    "table",
    "dict",
]

# Legitimately None / empty without related data -- just check access works.
NULLABLE_PROPERTIES = [
    "alias",
    "metadata",
    "elabs",
    "reaction",
    "scaffolds",
    "num_atoms_added",
]


@pytest.fixture
def compound(make_compound):
    """A registered ``Compound`` component (phenol)."""
    from designdb.components.compound import Compound

    return Compound(make_compound("c1ccccc1O"))


def test_not_null_properties(compound):
    for prop in NOT_NULL_PROPERTIES:
        assert getattr(compound, prop) is not None, f"{prop} is None"


def test_nullable_properties_do_not_raise(compound):
    for prop in NULLABLE_PROPERTIES:
        getattr(compound, prop)  # accessing must not raise


def test_core_values(compound):
    """Spot-check the computed chemistry for phenol (C6H5OH)."""
    assert compound.smiles == "c1ccccc1O"
    assert compound.inchikey == "ISWSIDIOOBJBQZ-UHFFFAOYSA-N"
    assert compound.name == compound.inchikey  # no alias -> falls back to inchikey
    assert compound.num_heavy_atoms == 7  # 6 C + 1 O
    assert compound.num_rings == 1
    assert compound.num_poses == 0
    assert compound.is_scaffold is False
