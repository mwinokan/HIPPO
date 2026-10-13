"""Subsite tests (SQLite tier).

Ported from the pre-refactor ``test_06_subsites.py`` (assign subsites from the
``CanonSites alias`` pose metadata field) and ``test_subsite.py`` (subsite
properties).
"""

import pytest

pytestmark = pytest.mark.sqlite


@pytest.fixture
def target(animal):
    """A dedicated target, so the subsites created here don't leak into the
    shared ``animal.target`` (``test_target.py`` expects it to have none)."""
    from designdb.models import TargetModel

    return TargetModel.objects.create(
        target_name=f"subsite-test-{TargetModel.objects.count()}",
        project=animal.target.project,
    )


@pytest.fixture
def poses(target, make_compound):
    """Three poses: two in subsite "sub-A", one in "sub-B", plus one without
    the metadata field."""
    from designdb.models import PoseModel
    from designdb.sets.pose import PoseSet

    compound = make_compound("c1ccc(Cl)cc1")  # chlorobenzene

    metadata = [
        {"CanonSites alias": "sub-A"},
        {"CanonSites alias": "sub-A"},
        {"CanonSites alias": "sub-B"},
        {"other": "field"},
    ]
    created = [
        PoseModel.objects.create(
            compound=compound,
            target=target,
            pose_alias=f"subsite-{i}",
            pose_metadata=m,
        )
        for i, m in enumerate(metadata)
    ]
    return PoseSet(PoseModel.objects.filter(pk__in=[p.pk for p in created]))


def test_set_subsites_from_metadata_field(target, poses):
    from designdb.models import SubsiteModel, SubsiteTagModel

    poses.set_subsites_from_metadata_field()

    subsites = SubsiteModel.objects.filter(
        target=target, subsite_name__in=["sub-A", "sub-B"]
    )
    assert set(subsites.values_list("subsite_name", flat=True)) == {"sub-A", "sub-B"}

    # the pose without the field is skipped, not failed
    tags = SubsiteTagModel.objects.filter(pose__in=poses.queryset)
    assert tags.count() == 3


def test_set_subsites_is_idempotent(poses):
    from designdb.models import SubsiteTagModel

    poses.set_subsites_from_metadata_field()
    poses.set_subsites_from_metadata_field()

    assert SubsiteTagModel.objects.filter(pose__in=poses.queryset).count() == 3


def test_subsite_properties(target, poses):
    from designdb.models import SubsiteModel

    poses.set_subsites_from_metadata_field()

    subsite = SubsiteModel.objects.get(target=target, subsite_name="sub-A")

    assert subsite.id is not None
    assert subsite.target == target
    assert subsite.subsite_name == "sub-A"
    assert subsite.posemodels.filter(pk__in=poses.queryset).count() == 2


def test_poseset_subsite_accessors(target, poses):
    from designdb.models import SubsiteModel

    poses.set_subsites_from_metadata_field()

    sub_a = SubsiteModel.objects.get(target=target, subsite_name="sub-A")

    assert len(poses.get_by_subsite(subsite=sub_a)) == 2
    assert len(poses(subsite=sub_a.pk)) == 2
    assert len(set(poses.subsite_ids)) == 2
    assert poses.num_subsites == 2

    df = poses.subsite_summary()
    assert df.loc[sub_a.pk, "num_poses"] == 2
