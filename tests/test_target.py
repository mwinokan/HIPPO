"""Target properties, modernized onto the SQLite ``animal`` conftest fixture.

Replaces the pre-refactor version (positional ``hippo.HIPPO('test', DB)``,
``animal.T1``, ``animal.db.close()``), which targeted the removed legacy API.
``animal.target`` is now a plain :class:`TargetModel`.
"""

import pytest

pytestmark = pytest.mark.sqlite


def test_target_identity(animal):
    """animal.target is the configured TargetModel, linked to its project."""
    target = animal.target

    assert target.pk is not None
    assert target.target_name == "test"
    # the project is created from the target_access_string (see conftest)
    assert target.project.project_name == "test-proposal"


def test_target_has_no_features_or_subsites_when_empty(animal):
    """A freshly-created target has no features/subsites until hits are loaded."""
    from designdb.models import FeatureModel, SubsiteModel

    target = animal.target

    assert FeatureModel.objects.filter(target=target).count() == 0
    assert SubsiteModel.objects.filter(target=target).count() == 0
