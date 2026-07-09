"""CompoundSet indexing tests (SQLite tier).

Regression for the same two ``__getitem__`` bugs fixed in PoseSet:
- a slice was passed straight to ``filter(pk__in=key)`` (errored);
- integer indexing didn't support negative indices (e.g. ``cset[-1]``).
"""

import pytest

pytestmark = pytest.mark.sqlite


@pytest.fixture
def compoundset(make_compound):
    """A CompoundSet of 10 distinct compounds."""
    from designdb.models import CompoundModel
    from designdb.sets.compound import CompoundSet

    smiles = ["C", "CC", "CCC", "CCCC", "CCCCC", "c1ccccc1", "CCO", "CCN", "CCCl", "CCBr"]
    ids = [make_compound(s).pk for s in smiles]
    return CompoundSet(CompoundModel.objects.filter(pk__in=ids))


def test_compoundset_slice_returns_positional_subset(compoundset):
    from designdb.sets.compound import CompoundSet

    sliced = compoundset[1:5]

    assert isinstance(sliced, CompoundSet)
    assert len(sliced) == 4
    assert list(sliced.ids) == list(compoundset.ids)[1:5]


def test_compoundset_slice_full_and_out_of_range(compoundset):
    assert len(compoundset[:]) == len(compoundset)
    assert len(compoundset[100:200]) == 0


def test_compoundset_slice_after_evaluation(compoundset):
    """Slicing must work even after the underlying queryset was evaluated.

    An evaluated queryset returns a *list of model instances* when sliced (from
    its result cache), which previously broke CompoundSet construction.
    """
    from designdb.sets.compound import CompoundSet

    list(compoundset)  # force-evaluate the underlying queryset (real-world usage)

    sliced = compoundset[1:5]

    assert isinstance(sliced, CompoundSet)
    assert len(sliced) == 4
    assert list(sliced.ids) == list(compoundset.ids)[1:5]


def test_compoundset_int_indexing_is_positional(compoundset):
    """cset[i] selects the i-th member by position, incl. negative indices."""
    from designdb.models import CompoundModel

    ids = list(compoundset.ids)

    first = compoundset[0]
    assert isinstance(first, CompoundModel)
    assert first.pk == ids[0]

    assert compoundset[1].pk == ids[1]
    assert compoundset[-1].pk == ids[-1]  # negative indexing -> last compound


def test_compoundset_int_index_out_of_range_raises(compoundset):
    with pytest.raises(IndexError):
        compoundset[len(compoundset)]
