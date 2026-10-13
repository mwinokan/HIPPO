"""Pose component tests (SQLite tier).

Ported from the pre-refactor ``test_pose.py`` property check. Properties that
need a protein structure on disk (``protein_system``, interactions) are left to
the integration tier, which has real Fragalysis data.
"""

import pytest

pytestmark = pytest.mark.sqlite


@pytest.fixture
def pose(animal, make_compound):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    from designdb.components.pose import Pose
    from designdb.models import PoseModel

    smiles = "c1ccc(N)cc1"  # aniline
    compound = make_compound(smiles)

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    mol = Chem.RemoveHs(mol)

    reference = PoseModel.objects.create(
        compound=compound, target=animal.target, pose_alias="pose-ref"
    )
    instance = PoseModel.objects.create(
        compound=compound,
        target=animal.target,
        pose_alias="pose-test",
        pose_smiles=smiles,
        pose_mol=mol,
        pose_metadata={"key": "value"},
        pose_reference=reference.pk,
    )
    return Pose(PoseModel.objects.get(pk=instance.pk))


def test_pose_identity(pose):
    assert pose.id is not None
    assert pose.pk == pose.id
    assert pose.instance.pk == pose.id
    assert str(pose) == f"P{pose.id}"


def test_pose_delegates_to_model(animal, pose):
    assert pose.pose_alias == "pose-test"
    assert pose.target == animal.target
    assert pose.compound.compound_smiles == "c1ccc(N)cc1"
    assert pose.pose_metadata == {"key": "value"}


def test_pose_mol_roundtrip(pose):
    from rdkit import Chem

    mol = pose.mol

    assert mol is not None
    assert mol.GetNumHeavyAtoms() == 7
    assert mol.GetNumConformers() == 1
    assert Chem.MolToSmiles(mol) == Chem.CanonSmiles("c1ccc(N)cc1")


def test_pose_features(pose):
    assert len(pose.features) > 0


def test_pose_reference(pose):
    from designdb.components.pose import Pose

    assert isinstance(pose.reference, Pose)
    assert pose.reference.id == pose.reference_id
    assert pose.reference.reference is None


def test_pose_fingerprint_flag(pose):
    from designdb.models import PoseModel

    assert not pose.has_fingerprint

    pose.set_has_fingerprint(True)

    assert pose.has_fingerprint
    assert PoseModel.objects.get(pk=pose.id).pose_fingerprint == 1


def test_pose_without_protein_link_has_no_protein_system(pose):
    assert pose.protein_link is None
    assert pose.protein_system is None
