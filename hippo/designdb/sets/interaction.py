"""Sets of protein-ligand interactions (ORM-backed).

An :class:`.InteractionSet` wraps a set of :class:`.InteractionModel` rows (via a
Django queryset). Construct it through :meth:`.PoseModel.interactions` /
:meth:`.PoseSet.interactions` or the factories here.

Interaction *detection* and duplicate *resolution* live in
:class:`.InteractionService` (``services/interaction.py``); this module is just the
read/aggregate surface over already-stored interactions.
"""

import mcol
import mrich
from designdb.models import InteractionModel
from django.db.models import Count

# `df` columns: ORM field (model / joined feature) -> output column name
_DF_COLUMNS = {
    'id': 'id',
    'feature_id': 'feature_id',
    'pose_id': 'pose_id',
    'feature__target_id': 'target_id',
    'interaction_type': 'type',
    'feature__feature_family': 'prot_family',
    'interaction_family': 'lig_family',
    'feature__feature_residue_name': 'residue_name',
    'feature__feature_residue_number': 'residue_number',
    'feature__feature_chain_name': 'chain_name',
    'interaction_distance': 'distance',
    'interaction_angle': 'angle',
    'interaction_energy': 'energy',
    'interaction_prot_coord': 'prot_coord',
    'interaction_lig_coord': 'lig_coord',
    'feature__feature_atom_name': 'prot_atoms',
    'interaction_atom_id': 'lig_atoms',
}


class InteractionSet:
    """A set of :class:`.InteractionModel` rows.

    .. attention::

            Not constructed directly -- use :meth:`.PoseModel.interactions` /
            :meth:`.PoseSet.interactions`, or the factory classmethods here.
    """

    def __init__(self, indices: list | None = None) -> None:
        """InteractionSet initialisation from a list of :class:`.InteractionModel` IDs"""
        indices = indices or []
        if not isinstance(indices, list):
            indices = list(indices)
        self._indices = sorted({int(i) for i in indices})
        self._df = None
        self._qs = InteractionModel.objects.filter(pk__in=self._indices)

    ### FACTORIES

    @classmethod
    def from_pose(cls, pose: 'PoseModel | PoseSet') -> 'InteractionSet':
        """Construct from one or more poses.

        :param pose: a :class:`.PoseSet` (has ``.ids``) or a single
            :class:`.Pose`/:class:`.PoseModel` (has ``.id``)
        """
        if hasattr(pose, 'ids'):
            qs = InteractionModel.objects.filter(pose_id__in=list(pose.ids))
        else:
            qs = InteractionModel.objects.filter(pose_id=pose.id)
        return cls(list(qs.values_list('id', flat=True)))

    @classmethod
    def all(cls) -> 'InteractionSet':
        """Construct an :class:`.InteractionSet` for every interaction in the table."""
        return cls(list(InteractionModel.objects.values_list('pk', flat=True)))

    @classmethod
    def from_residue(
        cls,
        residue_number: int,
        chain: str | None = None,
        target: 'TargetModel | int' = 1,
    ) -> 'InteractionSet':
        """Interactions formed with a given protein residue (and optionally chain).

        :param residue_number: the residue number
        :param chain: the chain name, or ``None`` for any chain
        :param target: the :class:`.TargetModel` or its ID (defaults to ``1``)
        """
        from designdb.models import TargetModel

        if isinstance(target, TargetModel):
            target = target.id

        qs = InteractionModel.objects.filter(
            feature__target_id=target,
            feature__feature_residue_number=residue_number,
        )
        if chain:
            qs = qs.filter(feature__feature_chain_name=chain)

        return cls(list(qs.values_list('id', flat=True)))

    ### PROPERTIES

    @property
    def queryset(self):
        """The underlying :class:`.InteractionModel` queryset"""
        return self._qs

    @property
    def indices(self) -> list[int]:
        """:class:`.InteractionModel` IDs in this set"""
        return self._indices

    @property
    def ids(self) -> list[int]:
        """:class:`.InteractionModel` IDs in this set"""
        return self._indices

    @property
    def types(self) -> list[str]:
        """Distinct interaction types in this set"""
        return list(self._qs.values_list('interaction_type', flat=True).distinct())

    @property
    def feature_ids(self) -> list[int]:
        """Distinct :class:`.FeatureModel` IDs interacted with"""
        return list(self._qs.values_list('feature_id', flat=True).distinct())

    @property
    def df(self) -> 'pandas.DataFrame':
        """DataFrame of the interactions, one row each."""
        if self._df is None:
            from pandas import DataFrame

            rows = list(self._qs.values(*_DF_COLUMNS))
            df = DataFrame(rows)
            if not df.empty:
                df = df.rename(columns=_DF_COLUMNS)
            self._df = df
        return self._df

    @property
    def _feature_counts(self) -> dict[int, int]:
        """Map of :class:`.FeatureModel` ID -> number of interactions with it"""
        return {
            row['feature']: row['n']
            for row in self._qs.values('feature').annotate(n=Count('id'))
        }

    @property
    def classic_fingerprint(self) -> dict:
        """Classic HIPPO fingerprint: :class:`.FeatureModel` ID -> interaction count."""
        return self.get_classic_fingerprint()

    @property
    def residue_number_chain_pairs(self) -> list[tuple]:
        """Distinct ``(residue_number, chain_name)`` pairs"""
        return list(
            self._qs.values_list(
                'feature__feature_residue_number', 'feature__feature_chain_name'
            ).distinct()
        )

    @property
    def type_residue_number_chain_triples(self) -> list[tuple]:
        """Distinct ``(interaction_type, residue_number, chain_name)`` triples"""
        return list(
            self._qs.values_list(
                'interaction_type',
                'feature__feature_residue_number',
                'feature__feature_chain_name',
            ).distinct()
        )

    @property
    def avg_num_residues_per_pose(self) -> float:
        """Mean number of distinct ``(residue, chain)`` contacts per pose"""
        from collections import defaultdict

        from numpy import mean

        d: dict[int, set] = defaultdict(set)
        for pose_id, res_num, chain in self._qs.values_list(
            'pose_id',
            'feature__feature_residue_number',
            'feature__feature_chain_name',
        ):
            d[pose_id].add((res_num, chain))
        return mean([len(v) for v in d.values()]) if d else 0

    @property
    def avg_num_interactions_per_pose(self) -> float:
        """Mean number of interactions per pose"""
        from collections import defaultdict

        from numpy import mean

        d: dict[int, int] = defaultdict(int)
        for (pose_id,) in self._qs.values_list('pose_id'):
            d[pose_id] += 1
        return mean(list(d.values())) if d else 0

    @property
    def avg_num_interaction_type_residue_pairs_per_pose(self) -> float:
        """Mean number of distinct ``(residue, type, chain)`` contacts per pose"""
        from collections import defaultdict

        from numpy import mean

        d: dict[int, set] = defaultdict(set)
        for pose_id, itype, res_num, chain in self._qs.values_list(
            'pose_id',
            'interaction_type',
            'feature__feature_residue_number',
            'feature__feature_chain_name',
        ):
            d[pose_id].add((res_num, itype, chain))
        return mean([len(v) for v in d.values()]) if d else 0

    @property
    def num_features(self) -> int:
        """Number of distinct protein :class:`.FeatureModel`\\ s interacted with"""
        return self._qs.values('feature').distinct().count()

    @property
    def avg_num_interactions_per_feature(self) -> float:
        """Mean number of interactions formed with each protein feature"""
        from numpy import mean

        counts = list(self._feature_counts.values())
        return mean(counts) if counts else 0

    @property
    def per_feature_count_hirsch(self) -> float:
        """h-index-like measure of how evenly features are interacted with"""
        from hirsch import hirsch

        counts = list(self._feature_counts.values())
        return hirsch(counts) if counts else 0

    ### METHODS

    def get_classic_fingerprint(self) -> dict:
        """Classic HIPPO fingerprint: :class:`.FeatureModel` ID -> interaction count."""
        return dict(self._feature_counts)

    def summary(self, families: bool = False) -> None:
        """Print a summary of this :class:`.InteractionSet`"""
        mrich.header(self)
        for i in self._qs.select_related('feature'):
            feature = i.feature
            s = (
                f'{i.interaction_type} '
                f'{feature.feature_residue_name}{feature.feature_residue_number}'
            )
            if families:
                s += f' {feature.feature_family} ~ {i.interaction_family}'
            mrich.var(s, f'{i.interaction_distance:.1f}', 'Å')

    ### DUNDERS

    def __len__(self) -> int:
        """The number of interactions in this set"""
        return len(self._indices)

    def __str__(self) -> str:
        """Unformatted command-line representation"""
        return f'{{I × {len(self)}}}'

    def __repr__(self) -> str:
        """ANSI formatted command-line representation"""
        return f'{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}'

    def __rich__(self) -> str:
        """Rich formatted command-line representation"""
        return f'[bold underline]{self}'

    def __iter__(self):
        """Iterate through the :class:`.InteractionModel` rows in this set"""
        return iter(self._qs)

    def __getitem__(self, key) -> 'InteractionModel | InteractionSet':
        """Index by position (int) or slice"""
        match key:
            case int():
                return InteractionModel.objects.get(pk=self._indices[key])
            case slice():
                return InteractionSet(self._indices[key])
            case _:
                raise NotImplementedError
