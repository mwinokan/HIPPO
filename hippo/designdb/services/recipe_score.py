"""Recipe scoring (user entry point: ``animal.scorers``).

A :class:`.Scorer` loads recipes from a directory of ``Recipe_*.json`` files and
evaluates them against weighted :class:`.Attribute` / :class:`.CustomAttribute`
objects: each attribute value is converted to a percentile (0-1) and combined by
weight. The score cache is written to ``{out_key}.json``.
"""

from pathlib import Path

import mrich
import numpy as np
import pandas as pd
from designdb.models import InteractionModel, PoseModel, ScaffoldModel
from designdb.recipe import Recipe, RecipeSet
from designdb.sets.compound import CompoundSet
from designdb.sets.interaction import InteractionSet
from designdb.sets.pose import PoseSet
from scipy.interpolate import interp1d

# columns of the internal score-cache DataFrame (besides the attribute columns)
DATA_COLUMNS = [
    'score',
    'price',
    'compound_ids',
    'pose_ids',
    'interaction_ids',
    'pose_metadata',
]


class Attribute:
    """A scoring attribute evaluated over the recipes of a :class:`.Scorer`.

    The raw value (``recipe.<key>``) is converted to a percentile in ``[0, 1]``
    across all recipes, optionally inverted, and scaled by ``weight``.
    """

    _type = 'Attribute'

    def __init__(
        self,
        scorer: 'Scorer',
        key: str,
        *,
        inverse: bool = False,
        weight: float = 1.0,
        bins: int = 100,
    ) -> None:
        self._scorer = scorer
        self._key = key
        self._inverse = inverse
        self._weight = weight
        self._bins = bins
        self._percentile_interpolator = None

    ### PROPERTIES

    @property
    def scorer(self) -> 'Scorer':
        return self._scorer

    @property
    def key(self) -> str:
        return self._key

    @property
    def inverse(self) -> bool:
        return self._inverse

    @property
    def bins(self) -> int:
        return self._bins

    @property
    def values(self) -> list[float]:
        """Attribute values across all recipes (computed/cached on access)."""
        col = self.scorer._data[self.key]
        null = col.isnull()
        if null.sum():
            for key in col[null].index.values:
                self.get_value(self.scorer.recipes[key], force=True)
            self.scorer._dump_json()
        return list(self.scorer._data[self.key].to_dict().values())

    @property
    def mean(self) -> float:
        return float(np.mean(self.values))

    @property
    def std(self) -> float:
        return float(np.std(self.values))

    @property
    def max(self) -> float:
        return max(self.values)

    @property
    def min(self) -> float:
        return min(self.values)

    @property
    def weight(self) -> float:
        return self._weight

    @weight.setter
    def weight(self, w) -> None:
        self.scorer._flag_weight_modification()
        self._weight = abs(w)
        self._inverse = self._inverse or (w < 0)

    @property
    def percentile_interpolator(self):
        """Interpolator mapping a value to its cumulative percentile."""
        if self._percentile_interpolator is None:
            count, bins_count = np.histogram(self.values, bins=self.bins)
            pdf = count / sum(count)
            cdf = np.cumsum(pdf)
            self._percentile_interpolator = interp1d(
                bins_count[1:], cdf, kind='linear', fill_value='extrapolate'
            )
        return self._percentile_interpolator

    ### METHODS

    def get_value(self, recipe: 'Recipe', force: bool = False) -> float:
        """Get (and cache) this attribute's raw value for ``recipe``."""
        if not force:
            cached = self.scorer._data[self.key][recipe.hash]
            if cached is not None:
                return cached
        value = getattr(recipe, self.key)
        self.scorer._data.at[recipe.hash, self.key] = value
        return value

    def unweighted(self, recipe: 'Recipe') -> float:
        """Percentile score (0-1) for ``recipe``."""
        value = self.get_value(recipe)
        score = float(self.percentile_interpolator(value))
        if self.inverse:
            score = 1 - score
        return score

    def __call__(self, recipe: 'Recipe') -> float:
        """Weighted score for ``recipe``."""
        if not self.weight:
            return 0.0
        return self.weight * self.unweighted(recipe)

    def __str__(self) -> str:
        return (
            f'{self._type}("{self.key}", weight={self.weight:.2f}, '
            f'inverse={self.inverse})'
        )

    def __repr__(self) -> str:
        import mcol

        return f'{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}'

    def __rich__(self) -> str:
        return f'[bold underline]{self}'


class CustomAttribute(Attribute):
    """A scoring attribute whose value is computed by a user-supplied function."""

    _type = 'CustomAttribute'

    def __init__(self, scorer: 'Scorer', key: str, function) -> None:
        self._function = function
        super().__init__(scorer=scorer, key=key)

    def get_value(self, recipe: 'Recipe', force: bool = False) -> float:
        if not force:
            cached = self.scorer._data[self.key][recipe.hash]
            if cached is not None:
                return cached
        value = self._function(recipe)
        self.scorer._data.at[recipe.hash, self.key] = value
        return value


class Scorer:
    """Score a set of recipes against weighted attributes."""

    def __init__(
        self,
        directory: 'str | Path',
        *,
        pattern: str = '*.json',
        attributes: list[str] | None = None,
        populate: bool = True,
        load_cache: bool = True,
        allowed_poses: 'PoseSet | list[int] | None' = None,
        out_key: str = 'scorer',
    ) -> None:
        """Scorer initialisation"""

        self._out_key = out_key

        if allowed_poses is None:
            self._allowed_pose_ids = None
        elif isinstance(allowed_poses, PoseSet):
            self._allowed_pose_ids = set(allowed_poses.ids)
        else:
            self._allowed_pose_ids = set(allowed_poses)

        self._recipes = RecipeSet(directory, pattern=pattern)

        self._attributes = {}
        for key in attributes or []:
            self._attributes[key] = Attribute(self, key)

        self._data = pd.DataFrame(
            index=self._recipes.keys(),
            columns=DATA_COLUMNS + self.attribute_keys,
        )
        self._data = self._data.replace({np.nan: None})

        if populate:
            if load_cache and self.json_path.exists():
                self._load_json()
            else:
                self._populate_query_cache()
            self._populate_recipe_child_sets()

        self.weights = 1.0

    ### FACTORIES

    @classmethod
    def default(
        cls,
        directory: 'str | Path',
        *,
        pattern: str = '*.json',
        skip: list[str] | None = None,
        load_cache: bool = True,
        subsites: bool = True,
        allowed_poses: 'PoseSet | list[int] | None' = None,
        out_key: str = 'scorer',
    ) -> 'Scorer':
        """Create a :class:`.Scorer` with the default attributes (minus ``skip``)."""

        self = cls.__new__(cls)

        standard = [k for k, v in DEFAULT_ATTRIBUTES.items() if v['type'] == 'standard']

        self.__init__(
            directory=directory,
            pattern=pattern,
            attributes=standard,
            populate=False,
            allowed_poses=allowed_poses,
            out_key=out_key,
        )

        skip = list(skip or [])

        # drop metrics whose data isn't present (ORM checks, no db)
        if not InteractionModel.objects.exists():
            mrich.warning('No interactions in DB, skipping related metrics')
            skip += ['interaction_count', 'interaction_balance']
        if not PoseModel.objects.exists():
            mrich.warning('No poses in DB, skipping related metrics')
            skip += [
                'num_inspirations',
                'num_inspiration_sets',
                'avg_energy_score',
                'avg_distance_score',
            ]
        if not ScaffoldModel.objects.exists():
            mrich.warning('No scaffold entries in DB, skipping related metrics')
            skip += ['num_scaffolds', 'num_scaffolds_elaborated', 'elaboration_balance']

        for key, attribute in [
            (k, v) for k, v in DEFAULT_ATTRIBUTES.items() if v['type'] == 'custom'
        ]:
            if key in skip:
                continue
            if not subsites and 'subsite' in key:
                continue
            self.add_custom_attribute(
                key, attribute['function'], weight_reset_warning=False
            )

        if load_cache and self.json_path.exists():
            self._load_json()
        else:
            self._populate_query_cache()
        self._populate_recipe_child_sets()

        # normalise weights from the DEFAULT_ATTRIBUTES config
        wsum = sum(abs(d['weight']) for d in DEFAULT_ATTRIBUTES.values())
        for attribute in self.attributes:
            attribute.weight = DEFAULT_ATTRIBUTES[attribute.key]['weight'] / wsum

        return self

    ### PROPERTIES

    @property
    def num_recipes(self) -> int:
        return len(self._recipes)

    @property
    def attributes(self) -> list:
        return list(self._attributes.values())

    @property
    def attribute_keys(self) -> list[str]:
        return list(self._attributes.keys())

    @property
    def recipes(self) -> 'RecipeSet':
        return self._recipes

    @property
    def num_attributes(self) -> int:
        return len(self._attributes)

    @property
    def weights(self) -> list[float]:
        return [a.weight for a in self.attributes]

    @weights.setter
    def weights(self, ws) -> None:
        self._flag_weight_modification()
        if isinstance(ws, (int, float)):
            ws = [ws] * self.num_attributes
        ws = list(ws)
        wsum = sum(abs(w) for w in ws) or 1.0
        for a, w in zip(self.attributes, ws, strict=False):
            a.weight = w / wsum

    @property
    def json_path(self) -> 'Path':
        """Path the score cache is written to (derived from ``out_key``)."""
        return Path(f'{self._out_key}.json')

    @property
    def score_dict(self) -> dict:
        """Scores keyed by recipe hash (computed/cached on access)."""
        col = self._data['score']
        null = col.isnull()
        if null.sum():
            mrich.debug('Calculating scores...')
            for key in col[null].index.values:
                self._data.at[key, 'score'] = self.score(self.recipes[key])
            self._dump_json()
        return self._data['score'].to_dict()

    @property
    def scores(self) -> list[float]:
        return list(self.score_dict.values())

    @property
    def best(self) -> 'Recipe':
        """Highest-scoring recipe."""
        return self.top(1)

    ### METHODS

    def add_custom_attribute(
        self, key: str, function, weight_reset_warning: bool = True
    ) -> 'CustomAttribute':
        """Add a custom scoring attribute computed by ``function(recipe)``."""
        ca = CustomAttribute(self, key, function)
        if key not in self._attributes:
            self._attributes[key] = ca
            if weight_reset_warning:
                mrich.warning('Attribute weights have been reset')
            self.weights = 1.0
            self._data[key] = None
        else:
            mrich.warning(f'Existing attribute with {key=}')
        return self._attributes[key]

    def score(self, recipe: 'Recipe', *, debug: bool = False) -> float:
        """Total weighted score for ``recipe``."""
        score = sum(attribute(recipe) for attribute in self.attributes)
        recipe._score = score
        return score

    def get_sorted_df(self) -> 'pd.DataFrame':
        """Score cache sorted by descending score."""
        _ = self.scores
        return self._data.sort_values(by='score', ascending=False)

    def top_keys(self, n: int) -> list[str]:
        return list(self.get_sorted_df().index[:n])

    def top(self, n: int) -> 'Recipe | list[Recipe]':
        """Top-``n`` scoring recipes (a single recipe when ``n == 1``)."""
        keys = self.top_keys(n)
        recipes = [self.recipes[key] for key in keys]
        return recipes[0] if n == 1 else recipes

    def plot(self, keys: list[str], budget: float | None = None):
        """Scatter plot of two attributes, coloured by score."""
        import plotly.express as px

        if len(keys) != 2:
            mrich.error('Only two keys supported')
            return None

        _ = self.scores

        df = self._data.drop(
            columns=['compound_ids', 'pose_ids', 'interaction_ids', 'pose_metadata']
        )
        df['score'] = pd.to_numeric(df['score'])

        for key in keys:
            if key not in df.columns:
                raise KeyError(f'no attribute/column named "{key}"')

        if budget:
            df = df[df['price'] < budget]

        df['hash'] = df.index.values
        return px.scatter(
            df, x=keys[0], y=keys[1], color='score', hover_data=list(df.columns)
        )

    ### INTERNALS

    def _flag_weight_modification(self) -> None:
        """Reset cached scores (weights changed)."""
        if hasattr(self, '_data'):
            self._data['score'] = None

    def _populate_query_cache(self) -> None:
        """Pre-fetch per-recipe compound/pose/interaction IDs + pose metadata (ORM)."""

        df = self._data

        # prices + combined compound IDs
        mrich.debug('Populating _data["compound_ids"]...')
        for recipe in self.recipes:
            df.at[recipe.hash, 'price'] = recipe.price.amount
            df.at[recipe.hash, 'compound_ids'] = recipe.combined_compound_ids

        # compound -> pose IDs
        all_compound_ids = set().union(*(set(ids) for ids in df['compound_ids'] if ids))
        mrich.debug(f'Getting poses for {len(all_compound_ids)} compounds')
        compound_pose_map: dict[int, set] = {}
        for c_id, p_id in PoseModel.objects.filter(
            compound_id__in=all_compound_ids
        ).values_list('compound_id', 'id'):
            compound_pose_map.setdefault(c_id, set()).add(p_id)

        mrich.debug('Populating _data["pose_ids"]...')
        for key in df.index.values:
            pose_ids = set()
            for c_id in df['compound_ids'][key]:
                ids = compound_pose_map.get(c_id, set())
                if self._allowed_pose_ids is not None:
                    ids = {i for i in ids if i in self._allowed_pose_ids}
                pose_ids |= ids
            df.at[key, 'pose_ids'] = pose_ids

        # pose -> interaction IDs
        all_pose_ids = set().union(*(set(ids) for ids in df['pose_ids'] if ids))
        mrich.debug(f'Getting interactions for {len(all_pose_ids)} poses')
        pose_interaction_map: dict[int, set] = {}
        if all_pose_ids:
            for p_id, i_id in InteractionModel.objects.filter(
                pose_id__in=all_pose_ids
            ).values_list('pose_id', 'id'):
                pose_interaction_map.setdefault(p_id, set()).add(i_id)

        mrich.debug('Populating _data["interaction_ids"]...')
        for key in df.index.values:
            interaction_ids = set()
            for p_id in df['pose_ids'][key]:
                interaction_ids |= pose_interaction_map.get(p_id, set())
            df.at[key, 'interaction_ids'] = interaction_ids

        # pose -> metadata
        mrich.debug(f'Getting metadata for {len(all_pose_ids)} poses')
        metadata_map: dict[int, dict] = {}
        if all_pose_ids:
            for p_id, meta in PoseModel.objects.filter(pk__in=all_pose_ids).values_list(
                'id', 'pose_metadata'
            ):
                metadata_map[p_id] = meta or {}

        mrich.debug('Populating _data["pose_metadata"]...')
        for key in df.index.values:
            df.at[key, 'pose_metadata'] = {
                p_id: metadata_map.get(p_id, {}) for p_id in df['pose_ids'][key]
            }

    def _populate_recipe_child_sets(self) -> None:
        """Inject the pre-fetched child sets onto each recipe's caches."""
        mrich.debug('Populating recipe caches')
        for key, recipe in self.recipes.items():
            row = self._data.loc[key]

            if recipe._combined_compounds is None:
                recipe._combined_compounds = CompoundSet(list(row['compound_ids']))
            if recipe._poses is None:
                recipe._poses = PoseSet(
                    PoseModel.objects.filter(pk__in=list(row['pose_ids']))
                )
            if recipe._interactions is None:
                recipe._interactions = InteractionSet(list(row['interaction_ids']))
            if recipe._poses._metadata_dict is None:
                recipe._poses._metadata_dict = row['pose_metadata']

    def _dump_json(self) -> None:
        path = self.json_path
        if path.parent and not path.parent.exists():
            path.parent.mkdir(parents=True)
        mrich.writing(path)
        self._data.to_json(path)

    def _load_json(self) -> None:
        path = self.json_path
        mrich.reading(path)
        cached = pd.read_json(path, orient='columns')

        if set(cached.columns) != set(self._data.columns):
            raise ValueError("Cached score JSON columns don't match expectation")
        if set(self._data.index.values) - set(cached.index.values):
            raise ValueError('Cached score JSON is missing recipes')

        self._data = cached.replace({np.nan: None})

    ### DUNDERS

    def __str__(self) -> str:
        return f'Scorer(#recipes={self.num_recipes})'

    def __repr__(self) -> str:
        import mcol

        return f'{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}'

    def __rich__(self) -> str:
        return f'[bold underline]{self}'


DEFAULT_ATTRIBUTES = {
    'num_scaffolds': dict(
        type='custom',
        weight=1.0,
        function=lambda r: r.combined_compounds.count_by_tag(tag='Syndirella scaffold'),
        description='Number of Syndirella scaffold compounds. Higher is better.',
    ),
    'num_compounds': dict(
        type='standard',
        weight=1.0,
        description='Number of product compounds. Higher is better.',
    ),
    'num_scaffolds_elaborated': dict(
        type='custom',
        weight=1.0,
        function=lambda r: r.combined_compounds.num_scaffolds_elaborated,
        description='Number of scaffolds with >=1 elaboration. Higher is better.',
    ),
    'elaboration_balance': dict(
        type='custom',
        weight=1.0,
        function=lambda r: r.combined_compounds.elaboration_balance,
        description='Evenness of scaffold elaboration (h-index). Higher is better.',
    ),
    'num_inspirations': dict(
        type='custom',
        weight=1.0,
        function=lambda r: r.poses.num_inspirations,
        description='Number of unique fragment inspirations. Higher is better.',
    ),
    'num_inspiration_sets': dict(
        type='custom',
        weight=1.0,
        function=lambda r: r.poses.num_inspiration_sets,
        description='Number of unique fragment combinations. Higher is better.',
    ),
    'interaction_count': dict(
        type='custom',
        weight=1.0,
        function=lambda r: r.interactions.num_features,
        description='Number of protein features interacted with. Higher is better.',
    ),
    'interaction_balance': dict(
        type='custom',
        weight=1.0,
        function=lambda r: r.interactions.per_feature_count_hirsch,
        description='Evenness of interactions across features. Higher is better.',
    ),
    'num_subsites': dict(
        type='custom',
        weight=1.0,
        function=lambda r: r.poses.num_subsites,
        description='Number of subsites occupied. Higher is better.',
    ),
    'subsite_balance': dict(
        type='custom',
        weight=1.0,
        function=lambda r: r.poses.subsite_balance,
        description='Evenness of subsite occupancy. Higher is better.',
    ),
    'avg_distance_score': dict(
        type='custom',
        weight=1.0,
        function=lambda r: r.poses.avg_distance_score,
        description='Average pose distance score.',
    ),
    'avg_energy_score': dict(
        type='custom',
        weight=1.0,
        function=lambda r: r.poses.avg_energy_score,
        description='Average pose energy score.',
    ),
}
