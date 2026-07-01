"""Client-side accessors (managers) for the HIPPO user-facing API.

Each manager is bound to a :class:`.HIPPO` instance and groups a family of
construction entry points (recipes, ingredients, routes, scorers, generators),
validating input and delegating to the backend services/sets. Together they form
the client surface of the eventual client/backend split (see CLAUDE.md).

Access them via the corresponding :class:`.HIPPO` properties: ``animal.recipes``,
``animal.ingredients``, ``animal.routes``, ``animal.scorers``, ``animal.generators``.
"""

from typing import TYPE_CHECKING

from .recipe import Recipe
from .services.generation import (
    RandomRecipeGenerator,
    RandomRecipeSelectionGenerator,
    RandomSelectionGenerator,
)
from .services.recipe import RecipeService
from .services.recipe_score import Scorer
from .sets.compound import CompoundSet
from .sets.ingredient import IngredientSet
from .sets.reaction import ReactionSet
from .sets.route import RouteSet

if TYPE_CHECKING:
    from .animal import HIPPO


class RecipeManager:
    """Build recipes from compounds/reactions/reactants (via ``animal.recipes``)."""

    def __init__(self, animal: 'HIPPO') -> None:
        self._animal = animal

    def from_compounds(self, compounds: CompoundSet, **kwargs):
        """Build recipe(s) to synthesise a :class:`.CompoundSet`.

        See :meth:`.RecipeService.from_compounds` for keyword arguments.
        """
        if not isinstance(compounds, CompoundSet):
            raise TypeError(f'compounds must be a CompoundSet, got {type(compounds)}')
        return RecipeService.from_compounds(compounds, **kwargs)

    def from_reactions(self, reactions: ReactionSet, **kwargs):
        """Build recipe(s) from a :class:`.ReactionSet`.

        See :meth:`.RecipeService.from_reactions` for keyword arguments.
        """
        if not isinstance(reactions, ReactionSet):
            raise TypeError(f'reactions must be a ReactionSet, got {type(reactions)}')
        return RecipeService.from_reactions(reactions, **kwargs)

    def from_reactants(self, reactants: 'CompoundSet | IngredientSet', **kwargs):
        """Build the maximal recipe reachable from a set of reactants.

        See :meth:`.RecipeService.from_reactants` for keyword arguments.
        """
        if not isinstance(reactants, (CompoundSet, IngredientSet)):
            raise TypeError(
                'reactants must be a CompoundSet or IngredientSet, '
                f'got {type(reactants)}'
            )
        return RecipeService.from_reactants(reactants, **kwargs)

    def from_json(self, path, **kwargs) -> 'Recipe':
        """Load a serialised :class:`.Recipe` from a JSON file.

        See :meth:`.Recipe.from_json` for keyword arguments (``data``,
        ``clear_quotes``, ``debug``).
        """
        return Recipe.from_json(path, **kwargs)


class IngredientManager:
    """Build :class:`.IngredientSet`\\ s (via ``animal.ingredients``)."""

    def __init__(self, animal: 'HIPPO') -> None:
        self._animal = animal

    def from_compounds(self, compounds: 'CompoundSet | None' = None, **kwargs):
        """Build an :class:`.IngredientSet` from a :class:`.CompoundSet` (or IDs).

        See :meth:`.IngredientSet.from_compounds` for keyword arguments.
        """
        if compounds is not None and not isinstance(compounds, CompoundSet):
            raise TypeError(f'compounds must be a CompoundSet, got {type(compounds)}')
        return IngredientSet.from_compounds(compounds=compounds, **kwargs)


class RouteManager:
    """Build :class:`.RouteSet`\\ s (via ``animal.routes``)."""

    def __init__(self, animal: 'HIPPO') -> None:
        self._animal = animal

    def from_product_ids(self, ids: 'CompoundSet | list[int]', *, progress=True):
        """Build a :class:`.RouteSet` of stored routes to the given products.

        :param ids: product :class:`.CompoundModel` IDs (or a :class:`.CompoundSet`)
        :param progress: show a progress bar while building
        """
        if isinstance(ids, CompoundSet):
            ids = ids.ids
        return RouteSet.from_product_ids(ids, progress=progress)


class ScorerManager:
    """Build recipe :class:`.Scorer`\\ s (via ``animal.scorers``)."""

    def __init__(self, animal: 'HIPPO') -> None:
        self._animal = animal

    def default(self, directory, **kwargs) -> Scorer:
        """Create a :class:`.Scorer` with the default attributes.

        See :meth:`.Scorer.default` for keyword arguments (``skip``,
        ``load_cache``, ``allowed_poses``, ``out_key``, ...).
        """
        return Scorer.default(directory, **kwargs)

    def create(self, directory, **kwargs) -> Scorer:
        """Create a :class:`.Scorer` with explicit attributes.

        See :class:`.Scorer` for keyword arguments.
        """
        return Scorer(directory, **kwargs)


class GeneratorManager:
    """Build random recipe/selection generators (via ``animal.generators``)."""

    def __init__(self, animal: 'HIPPO') -> None:
        self._animal = animal

    @staticmethod
    def _check(route_pool, compounds) -> None:
        if route_pool is not None and not isinstance(route_pool, RouteSet):
            raise TypeError(f'route_pool must be a RouteSet, got {type(route_pool)}')
        if compounds is not None and not isinstance(compounds, CompoundSet):
            raise TypeError(f'compounds must be a CompoundSet, got {type(compounds)}')

    def random_recipe(
        self,
        *,
        out_key: str,
        route_pool=None,
        suppliers=None,
        max_lead_time=None,
        start_with=None,
    ) -> RandomRecipeGenerator:
        """A generator that samples synthetic :class:`.Route`\\ s.

        :param out_key: base path/key for the generator's output files
        """
        self._check(route_pool, None)
        return RandomRecipeGenerator(
            out_key=out_key,
            route_pool=route_pool,
            suppliers=suppliers,
            max_lead_time=max_lead_time,
            start_with=start_with,
        )

    def random_selection(
        self,
        *,
        out_key: str,
        compounds=None,
        suppliers=None,
        amount: float = 1.0,
        max_lead_time=None,
        start_with=None,
    ) -> RandomSelectionGenerator:
        """A generator that samples (catalogue) compound selections.

        :param out_key: base path/key for the generator's output files
        """
        self._check(None, compounds)
        return RandomSelectionGenerator(
            out_key=out_key,
            compounds=compounds,
            suppliers=suppliers,
            amount=amount,
            max_lead_time=max_lead_time,
            start_with=start_with,
        )

    def random_recipe_selection(
        self,
        *,
        out_key: str,
        route_pool=None,
        compounds=None,
        suppliers=None,
        amount: float = 1.0,
        max_lead_time=None,
        start_with=None,
    ) -> RandomRecipeSelectionGenerator:
        """A generator combining routes and compound selections.

        :param out_key: base path/key for the generator's output files
        """
        self._check(route_pool, compounds)
        return RandomRecipeSelectionGenerator(
            out_key=out_key,
            route_pool=route_pool,
            compounds=compounds,
            suppliers=suppliers,
            amount=amount,
            max_lead_time=max_lead_time,
            start_with=start_with,
        )
