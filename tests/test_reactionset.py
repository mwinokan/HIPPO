"""ReactionSet indexing tests (SQLite tier).

Regression for the slice bug shared with Pose/CompoundSet: ``__getitem__`` for a
slice passed the ``slice`` object straight to ``filter(pk__in=key)`` (errored)
instead of positionally slicing the members.
"""

import pytest

pytestmark = pytest.mark.sqlite


@pytest.fixture
def reactionset(make_compound):
    """A ReactionSet of 10 reactions (one per distinct product compound)."""
    from designdb.models import ReactionModel
    from designdb.sets.reaction import ReactionSet

    smiles = ["C", "CC", "CCC", "CCCC", "CCCCC", "c1ccccc1", "CCO", "CCN", "CCCl", "CCBr"]
    reactions = [
        ReactionModel.objects.get_or_create(
            product_compound=make_compound(s), reaction_type="test"
        )[0]
        for s in smiles
    ]
    return ReactionSet(ReactionModel.objects.filter(pk__in=[r.pk for r in reactions]))


def test_reactionset_slice_returns_positional_subset(reactionset):
    from designdb.sets.reaction import ReactionSet

    sliced = reactionset[1:5]

    assert isinstance(sliced, ReactionSet)
    assert len(sliced) == 4
    assert list(sliced.ids) == list(reactionset.ids)[1:5]


def test_reactionset_slice_full_and_out_of_range(reactionset):
    assert len(reactionset[:]) == len(reactionset)
    assert len(reactionset[100:200]) == 0


def test_reactionset_slice_after_evaluation(reactionset):
    """Slicing must work even after the underlying queryset was evaluated."""
    from designdb.sets.reaction import ReactionSet

    list(reactionset)  # force-evaluate the underlying queryset

    sliced = reactionset[1:5]

    assert isinstance(sliced, ReactionSet)
    assert len(sliced) == 4
    assert list(sliced.ids) == list(reactionset.ids)[1:5]
