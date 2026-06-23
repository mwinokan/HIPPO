"""Pose component wrapping a :class:`.PoseModel`.

A :class:`.Pose` is a particular conformer of a :class:`.Compound` in a protein
environment. This component wraps the ORM :class:`.PoseModel`, exposing the ligand
molecule, the protein structure, and interaction-fingerprinting entry points.

Missing attributes are delegated to the wrapped :class:`.PoseModel`, so model
fields/relations (``pose_alias``, ``tags``, ``inspirations``, ``save`` …) remain
accessible.
"""

import mcol
import molparse as mp
from designdb.models import PoseModel


class Pose:
    """A conformer of a :class:`.Compound` within a protein environment."""

    def __init__(self, instance: 'PoseModel'):
        """Pose initialisation"""
        self._instance = instance
        self._protein_system = None

    def __getattr__(self, key: str):
        """Delegate unknown attributes to the wrapped :class:`.PoseModel`."""
        # guard internal attributes to avoid recursion before _instance is set
        if key.startswith('_'):
            raise AttributeError(key)
        return getattr(self._instance, key)

    ### PROPERTIES

    @property
    def instance(self) -> 'PoseModel':
        """The wrapped :class:`.PoseModel`"""
        return self._instance

    @property
    def id(self) -> int:
        """The pose's database ID"""
        return self._instance.id

    @property
    def pk(self) -> int:
        """The pose's primary key"""
        return self._instance.pk

    @property
    def mol(self):
        """The pose's ligand ``rdkit.Chem.Mol`` (stored in the DB)"""
        return self._instance.pose_mol

    @property
    def protein_link(self) -> str | None:
        """Path/link to the pose's protein structure (PDB)"""
        return self._instance.protein_link

    @property
    def reference_id(self) -> int | None:
        """ID of the pose's protein reference pose, if any"""
        return self._instance.pose_reference

    @property
    def reference(self) -> 'Pose | None':
        """The pose's protein reference (another :class:`.Pose`), if any"""
        ref_id = self._instance.pose_reference
        if ref_id is None:
            return None
        return Pose(PoseModel.objects.get(pk=ref_id))

    @property
    def protein_system(self) -> 'mp.System | None':
        """The pose's protein ``molparse.System`` (parsed from the PDB)"""
        if self._protein_system is None:
            link = self.protein_link
            if link and str(link).endswith('.pdb'):
                self._protein_system = mp.parse(link, verbosity=False).protein_system
        return self._protein_system

    @protein_system.setter
    def protein_system(self, system) -> None:
        """Set the pose's protein ``molparse.System``"""
        self._protein_system = system

    @property
    def features(self) -> list:
        """The pose ligand's ``molparse`` features"""
        return mp.rdkit.features_from_mol(self.mol)

    @property
    def has_fingerprint(self) -> bool:
        """Whether this pose has had its interactions fingerprinted"""
        return bool(self._instance.pose_fingerprint)

    ### METHODS

    def set_has_fingerprint(self, fp: bool, commit: bool = True) -> None:
        """Record whether this pose has been fingerprinted.

        :param fp: fingerprint state
        :param commit: persist to the database (Default value = True)
        """
        assert isinstance(fp, bool)
        self._instance.pose_fingerprint = int(fp)
        if commit:
            self._instance.save(update_fields=['pose_fingerprint'])

    def calculate_interactions(self, **kwargs) -> None:
        """Enumerate valid interactions between this pose's ligand and protein.

        Delegates to :meth:`.InteractionService.calculate`. See that method for
        keyword arguments.
        """
        from designdb.services.interaction import InteractionService

        return InteractionService.calculate(self, **kwargs)

    ### DUNDERS

    def __str__(self) -> str:
        """Plain string representation"""
        return f'P{self.id}'

    def __repr__(self) -> str:
        """ANSI formatted string representation"""
        return f'{mcol.bold}{mcol.underline}{str(self)}{mcol.unbold}{mcol.ununderline}'

    def __rich__(self) -> str:
        """Representation for mrich"""
        return f'[bold underline]{str(self)}'

    def __eq__(self, other) -> bool:
        """Equality by pose ID"""
        if isinstance(other, Pose):
            return self.id == other.id
        if isinstance(other, PoseModel):
            return self.id == other.id
        return NotImplemented

    def __hash__(self) -> int:
        """Hash by pose ID"""
        return hash(('Pose', self.id))
