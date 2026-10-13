"""Tag accessor tests (SQLite tier).

Ported from the pre-refactor ``test_tags.py``, which checked the unique tags and
the tag summary of the old ``animal.tags`` table. Tags now live on the sets:
``PoseSet.tags`` / ``CompoundSet.tags`` / ``CompoundSet.tag_summary``. Tag
*resolution* during ingestion is covered by ``test_tag_resolution.py``.
"""

import pytest

pytestmark = pytest.mark.sqlite


@pytest.fixture
def compounds(make_compound):
    from designdb.models import CompoundModel
    from designdb.sets.compound import CompoundSet

    created = [make_compound(s) for s in ("CCN", "CCCN", "CCCCN")]
    return CompoundSet(CompoundModel.objects.filter(pk__in=[c.pk for c in created]))


@pytest.fixture
def poses(animal, compounds):
    from designdb.models import PoseModel
    from designdb.sets.pose import PoseSet

    created = [
        PoseModel.objects.create(
            compound_id=cid, target=animal.target, pose_alias=f"tags-{cid}"
        )
        for cid in compounds.ids
    ]
    return PoseSet(PoseModel.objects.filter(pk__in=[p.pk for p in created]))


def test_compoundset_tags(animal, compounds):
    compounds.add_tag("tags-all", target=animal.target)
    compounds[0:1].add_tag("tags-one", target=animal.target)

    assert compounds.tags == {"tags-all", "tags-one"}
    assert compounds.get_tags(target=animal.target) == {"tags-all", "tags-one"}


def test_compoundset_tag_summary(animal, compounds):
    compounds.add_tag("tags-summary-all", target=animal.target)
    compounds[0:1].add_tag("tags-summary-one", target=animal.target)

    df = compounds.tag_summary()

    assert df.loc["tags-summary-all", "num_compounds"] == 3
    assert df.loc["tags-summary-one", "num_compounds"] == 1


def test_poseset_tags(poses):
    poses.add_tag("tags-pose-all")
    poses[0:1].add_tag("tags-pose-one")

    assert set(poses.tags) == {"tags-pose-all", "tags-pose-one"}
    assert len(poses(tag="tags-pose-all")) == 3
    assert len(poses.get_by_tag("tags-pose-one")) == 1
