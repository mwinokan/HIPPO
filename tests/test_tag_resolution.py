"""Tag registration and resolution tests (SQLite tier).

The tag vocabulary (``compound_tags`` / ``pose_tags``) only grows through
``animal.register_compound_tag`` / ``register_pose_tag``, which need DesignDB
admin rights. Ingestion and the ``add_tag`` helpers *resolve* tag names against
existing rows and never create them.
"""

import pytest

pytestmark = pytest.mark.sqlite


def queryset_for(compound):
    """A one-member queryset, for wrapping a compound in a CompoundSet."""
    from designdb.models import CompoundModel

    return CompoundModel.objects.filter(pk=compound.pk)


def compound_tags_of(compound, target):
    from designdb.models import CompoundTagJunctionModel

    return set(
        CompoundTagJunctionModel.objects.filter(
            compound=compound, target=target
        ).values_list("compound_tag__compound_tag_name", flat=True)
    )


@pytest.fixture
def compound_tag(animal):
    """An existing compound-tag vocabulary row."""
    from designdb.models import CompoundTagModel

    tag, _ = CompoundTagModel.objects.get_or_create(compound_tag_name="res-known-cmpd")
    return tag


@pytest.fixture
def pose_tag(animal):
    """An existing pose-tag vocabulary row."""
    from designdb.models import PoseTagModel

    tag, _ = PoseTagModel.objects.get_or_create(pose_tag_name="res-known-pose")
    return tag


### resolve_tags


def test_resolve_tags_returns_existing_rows(compound_tag, pose_tag):
    from designdb.services.compound import CompoundTagService
    from designdb.services.pose import PoseTagService

    assert CompoundTagService.resolve_tags(["res-known-cmpd"]) == [compound_tag]
    assert PoseTagService.resolve_tags(["res-known-pose"]) == [pose_tag]


def test_resolve_tags_returns_a_list_not_a_lazy_queryset(compound_tag):
    """ingest_sdf iterates the result once per chunk; a queryset would re-query."""
    from designdb.services.compound import CompoundTagService

    assert isinstance(CompoundTagService.resolve_tags(["res-known-cmpd"]), list)


def test_resolve_tags_strips_whitespace(compound_tag):
    """Regression: rows were created stripped but selected unstripped."""
    from designdb.services.compound import CompoundTagService

    assert CompoundTagService.resolve_tags([" res-known-cmpd "]) == [compound_tag]


def test_resolve_tags_ignores_blanks_and_duplicates(compound_tag):
    from designdb.services.compound import CompoundTagService

    assert CompoundTagService.resolve_tags([]) == []
    assert CompoundTagService.resolve_tags(["", "  ", None]) == []
    assert CompoundTagService.resolve_tags(
        ["res-known-cmpd", "res-known-cmpd"]
    ) == [compound_tag]


def test_resolve_tags_rejects_none():
    from designdb.services.compound import CompoundTagService

    with pytest.raises(ValueError):
        CompoundTagService.resolve_tags(None)


def test_resolve_tags_raises_naming_every_missing_tag(compound_tag):
    from designdb.services.compound import CompoundTagService
    from designdb.utils import MissingTagError

    with pytest.raises(MissingTagError) as exc:
        CompoundTagService.resolve_tags(
            ["res-known-cmpd", "res-absent-a", "res-absent-b"]
        )

    message = str(exc.value)
    assert "res-absent-a" in message
    assert "res-absent-b" in message
    # the tag that does exist is not part of the complaint
    assert "res-known-cmpd" not in message


def test_resolve_tags_creates_nothing(animal):
    from designdb.models import CompoundTagModel
    from designdb.services.compound import CompoundTagService
    from designdb.utils import MissingTagError

    before = CompoundTagModel.objects.count()

    with pytest.raises(MissingTagError):
        CompoundTagService.resolve_tags(["res-never-created"])

    assert CompoundTagModel.objects.count() == before
    assert not CompoundTagModel.objects.filter(
        compound_tag_name="res-never-created"
    ).exists()


def test_resolve_tags_is_a_value_error():
    """MissingTagError subclasses ValueError, matching MethodService."""
    from designdb.utils import MissingTagError

    assert issubclass(MissingTagError, ValueError)


### ingestion


def test_load_sdf_rejects_unknown_tag_before_touching_the_file(animal, tmp_path):
    """Tags resolve before any file IO, so a bad tag beats a missing SDF."""
    from designdb.models import CompoundModel, PoseModel
    from designdb.utils import MissingTagError

    compounds_before = CompoundModel.objects.count()
    poses_before = PoseModel.objects.count()

    with pytest.raises(MissingTagError):
        animal.load_sdf(
            path=tmp_path / "does-not-exist.sdf",
            compound_tags=["res-absent-ingest"],
        )

    assert CompoundModel.objects.count() == compounds_before
    assert PoseModel.objects.count() == poses_before


### register_tag


def test_register_compound_tag_creates_then_retrieves(animal):
    from designdb.models import CompoundTagModel

    tag, created = animal.register_compound_tag(
        "reg-cmpd", description="a description", note="a note"
    )
    assert created
    assert (tag.compound_tag_description, tag.compound_tag_note) == (
        "a description",
        "a note",
    )

    again, created = animal.register_compound_tag("reg-cmpd", description="ignored")
    assert not created
    assert again == tag
    assert again.compound_tag_description == "a description"
    assert CompoundTagModel.objects.filter(compound_tag_name="reg-cmpd").count() == 1
    assert tag in animal.compound_tags


def test_register_pose_tag_creates_then_retrieves(animal):
    tag, created = animal.register_pose_tag("reg-pose")
    assert created

    again, created = animal.register_pose_tag("reg-pose")
    assert not created
    assert again == tag
    assert tag in animal.pose_tags


def test_register_tag_strips_whitespace(animal):
    """Registered names must match what resolve_tags looks up (stripped)."""
    from designdb.services.compound import CompoundTagService

    tag, _ = animal.register_compound_tag("  reg-stripped  ")
    assert tag.compound_tag_name == "reg-stripped"
    assert CompoundTagService.resolve_tags(["reg-stripped"]) == [tag]


@pytest.mark.parametrize("name", ["", "   ", None])
def test_register_tag_rejects_blank_names(animal, name):
    with pytest.raises(ValueError):
        animal.register_compound_tag(name)
    with pytest.raises(ValueError):
        animal.register_pose_tag(name)


def test_registered_tag_is_accepted_by_ingestion(animal):
    """The error for an unknown tag points at register_*_tag; doing so fixes it."""
    from designdb.services.compound import CompoundTagService
    from designdb.utils import MissingTagError

    with pytest.raises(MissingTagError, match="register_compound_tag"):
        CompoundTagService.resolve_tags(["reg-later"])

    tag, _ = animal.register_compound_tag("reg-later")

    assert CompoundTagService.resolve_tags(["reg-later"]) == [tag]


### permission refusal


def permission_denied():
    """A ProgrammingError as Django raises it for Postgres SQLSTATE 42501."""
    from django.db import ProgrammingError

    class InsufficientPrivilege(Exception):
        sqlstate = "42501"

    exc = ProgrammingError("permission denied for table compound_tags")
    exc.__cause__ = InsufficientPrivilege()
    return exc


def test_admin_write_explains_a_permission_refusal(animal):
    from designdb.utils import InsufficientPrivilegeError, admin_write

    with pytest.raises(InsufficientPrivilegeError) as exc_info:
        with admin_write("Registering compound tag 'x'"):
            raise permission_denied()

    message = str(exc_info.value)
    assert "Registering compound tag 'x' needs DesignDB admin rights" in message
    assert "DesignDB administrator" in message
    assert isinstance(exc_info.value, PermissionError)


def test_admin_write_passes_other_errors_through(animal):
    from django.db import ProgrammingError

    from designdb.utils import admin_write

    with pytest.raises(ProgrammingError, match="syntax error"):
        with admin_write("Registering compound tag 'x'"):
            raise ProgrammingError("syntax error at or near")


def test_register_tag_reports_a_permission_refusal(animal):
    from unittest.mock import patch

    from designdb.models import CompoundTagModel
    from designdb.utils import InsufficientPrivilegeError

    with patch.object(
        CompoundTagModel.objects, "get_or_create", side_effect=permission_denied()
    ):
        with pytest.raises(InsufficientPrivilegeError, match="reg-refused"):
            animal.register_compound_tag("reg-refused")


### add_tag


def test_compoundset_add_tag_uses_a_registered_tag(animal, make_compound, compound_tag):
    from designdb.sets.compound import CompoundSet

    compound = make_compound("c1ccc(Cl)cc1")  # chlorobenzene
    cset = CompoundSet(queryset_for(compound))

    cset.add_tag(compound_tag.compound_tag_name, target=animal.target)

    assert compound_tag.compound_tag_name in compound_tags_of(compound, animal.target)


def test_compoundset_add_tag_rejects_an_unknown_tag(animal, make_compound):
    from designdb.models import CompoundTagModel
    from designdb.sets.compound import CompoundSet
    from designdb.utils import MissingTagError

    compound = make_compound("c1ccc(Br)cc1")  # bromobenzene
    cset = CompoundSet(queryset_for(compound))

    with pytest.raises(MissingTagError, match="register_compound_tag"):
        cset.add_tag("res-unknown-set", target=animal.target)

    assert not CompoundTagModel.objects.filter(
        compound_tag_name="res-unknown-set"
    ).exists()


def test_compound_add_tag_rejects_an_unknown_tag(animal, make_compound):
    from designdb.components.compound import Compound
    from designdb.models import CompoundTagModel
    from designdb.utils import MissingTagError

    compound = Compound(make_compound("c1ccc(I)cc1"))  # iodobenzene

    with pytest.raises(MissingTagError, match="register_compound_tag"):
        compound.add_tag("res-unknown-component", target=animal.target)

    assert not CompoundTagModel.objects.filter(
        compound_tag_name="res-unknown-component"
    ).exists()


def make_pose_set(animal, make_compound, smiles, alias):
    from designdb.models import PoseModel
    from designdb.sets.pose import PoseSet

    pose = PoseModel.objects.create(
        compound=make_compound(smiles), target=animal.target, pose_alias=alias
    )
    return pose, PoseSet(PoseModel.objects.filter(pk=pose.pk))


def test_poseset_add_tag_accepts_an_existing_tag(animal, make_compound, pose_tag):
    """Regression: a bare PoseTagModel(...).save() hit the uc_pose_tag constraint."""
    pose, pset = make_pose_set(animal, make_compound, "c1ccc(C)cc1", "res-tag-pose")

    pset.add_tag(pose_tag.pose_tag_name)  # must not raise IntegrityError

    assert pose_tag.pose_tag_name in {t.pose_tag_name for t in pose.tags.all()}


def test_poseset_add_tag_rejects_an_unknown_tag(animal, make_compound):
    from designdb.models import PoseTagModel
    from designdb.utils import MissingTagError

    _, pset = make_pose_set(animal, make_compound, "c1ccc(CC)cc1", "res-unknown-pose")

    with pytest.raises(MissingTagError, match="register_pose_tag"):
        pset.add_tag("res-unknown-pose-tag")

    assert not PoseTagModel.objects.filter(
        pose_tag_name="res-unknown-pose-tag"
    ).exists()
