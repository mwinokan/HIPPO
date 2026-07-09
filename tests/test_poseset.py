"""PoseSet tests (SQLite tier).

Regression: ``PoseSet.__getitem__`` for a slice passed the ``slice`` object
straight to ``filter(pk__in=key)`` instead of positionally slicing the members
(e.g. ``poseset[1:10]`` errored).
"""

import pytest

pytestmark = pytest.mark.sqlite


@pytest.fixture
def poseset(animal, make_compound):
    """A PoseSet of 20 minimal poses.

    Uses a compound dedicated to this test so it doesn't perturb the pose counts
    of compounds used by other tests sharing the session database.
    """
    from designdb.models import PoseModel
    from designdb.sets.pose import PoseSet

    compound = make_compound("c1ccc(F)cc1")  # fluorobenzene
    target = animal.target

    poses = [
        PoseModel.objects.create(
            compound=compound, target=target, pose_alias=f"slice-{i}"
        )
        for i in range(20)
    ]
    return PoseSet(PoseModel.objects.filter(pk__in=[p.pk for p in poses]))


def test_poseset_slice_returns_positional_subset(poseset):
    from designdb.sets.pose import PoseSet

    sliced = poseset[1:10]

    assert isinstance(sliced, PoseSet)
    assert len(sliced) == 9
    # a slice selects members by position, not by pk
    assert list(sliced.ids) == list(poseset.ids)[1:10]


def test_poseset_slice_full_and_out_of_range(poseset):
    assert len(poseset[:]) == 20
    assert len(poseset[100:200]) == 0


def test_poseset_slice_after_evaluation(poseset):
    """Slicing must work even after the underlying queryset was evaluated."""
    from designdb.sets.pose import PoseSet

    list(poseset)  # force-evaluate the underlying queryset

    sliced = poseset[1:10]

    assert isinstance(sliced, PoseSet)
    assert len(sliced) == 9
    assert list(sliced.ids) == list(poseset.ids)[1:10]


def test_poseset_int_indexing_is_positional(poseset):
    """pset[i] selects the i-th member by position (not by pk)."""
    from designdb.components.pose import Pose

    ids = list(poseset.ids)

    first = poseset[0]
    assert isinstance(first, Pose)
    assert first.id == ids[0]

    assert poseset[1].id == ids[1]  # second pose, by position
    assert poseset[-1].id == ids[-1]  # negative indexing -> last pose


def test_poseset_int_index_out_of_range_raises(poseset):
    with pytest.raises(IndexError):
        poseset[len(poseset)]
