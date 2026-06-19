"""Service layer: random recipe / selection generators.

Ported from the legacy ``rgen`` module and modernised:

* **No legacy ``Database`` coupling.** The modern :class:`.Recipe` / :class:`.Route`
  / :class:`.IngredientSet` / :class:`.Price` are self-contained (ORM-backed), so
  the generators take pools + config only -- no ``db`` / ``db.path``. As there is no
  database path to derive output filenames from, ``out_key`` is now **required**.
* **File output is preserved.** Each :meth:`generate` writes the generated recipe to
  ``{recipe_dir}/Recipe_<hash>.json`` (via :meth:`.Recipe.write_json`), and each
  generator dumps its state to ``{out_key}_<gen>.json`` on construction -- matching
  legacy behaviour. ``generate`` also returns the :class:`.Recipe` so callers may
  collect them in memory.
* The three legacy generators share one add-within-budget loop here
  (:func:`_generate_recipe`) rather than duplicating it.

Layering: ``services -> components/sets``. User entry point: ``animal.generators``.
"""

import json
from pathlib import Path

import mrich
from designdb.components.price import Price
from designdb.components.recipe import Recipe, Route
from designdb.sets.compound import IngredientSet
from designdb.sets.route import RouteSet
from designdb.utils import dt_hash

# sentinel for "no limit" on a given count
_UNLIMITED = 10**9


def _generate_recipe(
    starting_recipe: 'Recipe',
    pool: list,
    *,
    budget: 'Price',
    suppliers: 'list | None',
    max_products: int,
    max_reactions: int,
    max_compounds: int,
    max_iter: int | None,
    shuffle: bool,
    debug: bool,
) -> 'tuple[Recipe, dict]':
    """Randomly add routes/compounds from ``pool`` to a recipe within ``budget``.

    ``pool`` items are :class:`.Route` objects (added via ``recipe += route``) or
    :class:`.Ingredient` objects (added to ``recipe.compounds``). Candidates whose
    product/compound is already present are skipped; an addition that exceeds the
    budget is reverted. Stops on budget/limit/pool-depletion or ``max_iter``.

    :returns: ``(recipe, stats)`` where ``stats`` has ``stop_reason``,
        ``iterations`` and the resolved ``max_iter``.
    """

    from random import shuffle as shuffle_func

    recipe = starting_recipe.copy()

    if suppliers is not None:
        recipe.reactants._supplier = suppliers
        recipe.compounds._supplier = suppliers

    pool = list(pool)
    if not pool:
        raise ValueError('Route/compound pool is empty!')

    if max_iter is None:
        max_iter = max_products + max_reactions
    max_iter = min(max_iter, len(pool))

    if shuffle:
        mrich.debug('Shuffling pool')
        shuffle_func(pool)

    old_recipe = recipe.copy()
    stop_reason = 'max iterations reached'
    iterations = 0

    for i in mrich.track(range(max_iter), prefix='Generating recipe...'):
        iterations = i
        candidate = pool.pop()

        candidate_compound = (
            candidate.product if isinstance(candidate, Route) else candidate
        )

        if (
            candidate_compound in recipe.products
            or candidate_compound in recipe.compounds
        ):
            continue

        if isinstance(candidate, Route):
            recipe += candidate
        else:
            recipe.compounds.add(candidate)

        new_price = recipe.price

        if not pool:
            stop_reason = 'Route/compound pool depleted'
            break

        if new_price > budget:
            recipe = old_recipe.copy()
            continue

        if len(recipe.reactions) > max_reactions:
            stop_reason = 'Max #reactions exceeded'
            break

        if len(recipe.products) > max_products:
            stop_reason = 'Max #products exceeded'
            break

        if len(recipe.compounds) > max_compounds:
            stop_reason = 'Max #compounds exceeded'
            break

        old_recipe = recipe.copy()

    mrich.success(f'{stop_reason}!')
    mrich.success(f'Completed after {iterations} iterations!')

    stats = {'stop_reason': stop_reason, 'iterations': iterations, 'max_iter': max_iter}
    return recipe, stats


class _GeneratorBase:
    """Shared state, I/O setup and representation for the random generators."""

    _suppliers: 'list | None' = None
    _max_lead_time: 'float | None' = None
    _starting_recipe: 'Recipe | None' = None
    _out_key: str | None = None
    _data_path: 'Path | None' = None
    _recipe_dir: 'Path | None' = None

    def _setup_io(
        self,
        out_key: str | None,
        data_suffix: str,
        dir_suffix: str,
        skip_directory_creation: bool,
    ) -> None:
        """Resolve output paths from ``out_key`` and create directories.

        :param out_key: base path/key for output files (required -- there is no
            database path to derive it from)
        :param data_suffix: suffix for the state JSON, e.g. ``'_rsgen.json'``
        :param dir_suffix: suffix for the recipe directory, e.g.
            ``'_recipes_and_selections'``
        :param skip_directory_creation: don't create a recipe directory (used for
            the inner generators of a composed generator)
        """
        if not out_key:
            raise ValueError(
                'out_key is required (used for output filenames; there is no '
                'database path to derive it from)'
            )

        self._out_key = out_key

        parent_dir = Path(out_key).parent
        if not parent_dir.exists():
            parent_dir.mkdir(parents=True)

        self._data_path = Path(f'{out_key}{data_suffix}')
        if self._data_path.exists():
            mrich.warning(f'Will overwrite existing data file: {self._data_path}')

        if skip_directory_creation:
            self._recipe_dir = None
        else:
            path = Path(f'{out_key}{dir_suffix}')
            if not path.exists():
                mrich.writing(f'{path}/')
                path.mkdir()
            self._recipe_dir = path

    def _dump_data(self, data: dict) -> None:
        """Write the generator state dict to ``data_path``."""
        if self._recipe_dir is not None:
            data['recipe_dir'] = str(self._recipe_dir.resolve())
        else:
            data['recipe_dir'] = None
        data['suppliers'] = self._suppliers
        data['starting_recipe'] = self._starting_recipe.get_dict(serialise_price=True)
        mrich.writing(self._data_path)
        json.dump(data, open(self._data_path, 'wt'), indent=4)

    def _write_recipe(self, recipe: 'Recipe', budget: 'Price', stats: dict, params: dict):
        """Write a generated recipe to ``{recipe_dir}/Recipe_<hash>.json``."""
        out_file = self._recipe_dir / f'Recipe_{dt_hash()}.json'
        metadict = {
            'gen_data_path': str(self._data_path.resolve()),
            'gen_recipe_dir': str(self._recipe_dir.resolve()),
            'gen_max_lead_time': self._max_lead_time,
            'gen_suppliers': self._suppliers,
            'gen_budget': budget.amount,
            'gen_currency': budget.currency,
            'gen_stop_reason': stats['stop_reason'],
            'gen_iterations': stats['iterations'],
            'gen_max_iter': stats['max_iter'],
            'gen_recipe_path': str(out_file.resolve()),
            **params,
        }
        recipe.write_json(out_file, extra=metadict)
        return out_file

    @property
    def suppliers(self) -> 'list | None':
        """Restrict quoting to these suppliers"""
        return self._suppliers

    @property
    def max_lead_time(self) -> 'float | None':
        """Maximum lead-time constraint"""
        return self._max_lead_time

    @property
    def starting_recipe(self) -> 'Recipe':
        """The recipe every generation starts from (a copy)"""
        return self._starting_recipe

    @property
    def out_key(self) -> str | None:
        """Base key for output filenames"""
        return self._out_key

    @property
    def recipe_dir(self) -> 'Path | None':
        """Directory generated recipe JSONs are written to"""
        return self._recipe_dir

    def __call__(self, *args, **kwargs) -> 'Recipe':
        """Generate a recipe (alias for :meth:`generate`)"""
        return self.generate(*args, **kwargs)

    def __repr__(self) -> str:
        return f'{type(self).__name__}(out_key={self._out_key!r})'


class RandomRecipeGenerator(_GeneratorBase):
    """Generate random recipes by sampling synthetic :class:`.Route`\\ s."""

    def __init__(
        self,
        *,
        out_key: str,
        suppliers: 'list | None' = None,
        route_pool: 'RouteSet | None' = None,
        max_lead_time: 'float | None' = None,
        start_with: 'Recipe | None' = None,
        skip_directory_creation: bool = False,
    ) -> None:
        mrich.debug('RandomRecipeGenerator.__init__()')
        self._suppliers = suppliers
        self._max_lead_time = max_lead_time
        self._starting_recipe = start_with or Recipe()
        self._route_pool = route_pool if route_pool is not None else RouteSet()

        self._setup_io(out_key, '_rgen.json', '_recipes', skip_directory_creation)
        self._dump_data({'route_pool': self._route_pool.get_dict()})

    @property
    def route_pool(self) -> 'RouteSet':
        """Pool of routes to sample from"""
        return self._route_pool

    def generate(
        self,
        budget: float = 10000,
        currency: str = 'EUR',
        *,
        max_products: int = 1000,
        max_reactions: int = 1000,
        max_iter: int | None = None,
        shuffle: bool = True,
        balance_clusters: bool = False,
        permitted_clusters=None,
        debug: bool = False,
    ) -> 'Recipe':
        """Generate a random recipe of routes within ``budget`` (also written to disk)."""
        if balance_clusters:
            raise NotImplementedError(
                'balance_clusters requires route clustering, which is not yet ported'
            )
        budget = Price(budget, currency)
        recipe, stats = _generate_recipe(
            self._starting_recipe,
            list(self._route_pool),
            budget=budget,
            suppliers=self._suppliers,
            max_products=max_products,
            max_reactions=max_reactions,
            max_compounds=_UNLIMITED,
            max_iter=max_iter,
            shuffle=shuffle,
            debug=debug,
        )
        self._write_recipe(
            recipe,
            budget,
            stats,
            {'gen_max_products': max_products, 'gen_max_reactions': max_reactions},
        )
        return recipe


class RandomSelectionGenerator(_GeneratorBase):
    """Generate random selections of (catalogue) compounds."""

    def __init__(
        self,
        *,
        out_key: str,
        suppliers: 'list | None' = None,
        compounds=None,
        amount: float = 1.0,
        max_lead_time: 'float | None' = None,
        start_with: 'Recipe | None' = None,
        skip_directory_creation: bool = False,
    ) -> None:
        mrich.debug('RandomSelectionGenerator.__init__()')
        self._suppliers = suppliers
        self._max_lead_time = max_lead_time
        self._amount = amount
        self._starting_recipe = start_with or Recipe()
        if compounds is None:
            self._compound_pool = IngredientSet()
        else:
            self._compound_pool = IngredientSet.from_compounds(
                compounds=compounds, amount=amount
            )

        self._setup_io(out_key, '_sgen.json', '_selections', skip_directory_creation)
        self._dump_data(
            {'amount': amount, 'compound_pool': self._compound_pool.get_dict()}
        )

    @property
    def compound_pool(self) -> 'IngredientSet':
        """Pool of compound ingredients to sample from"""
        return self._compound_pool

    def generate(
        self,
        budget: float = 10000,
        currency: str = 'EUR',
        *,
        max_compounds: int = 1000,
        max_iter: int | None = None,
        shuffle: bool = True,
        debug: bool = False,
    ) -> 'Recipe':
        """Generate a random compound selection within ``budget`` (also written to disk)."""
        if max_iter is None:
            max_iter = max_compounds * 3
        budget = Price(budget, currency)
        recipe, stats = _generate_recipe(
            self._starting_recipe,
            list(self._compound_pool),
            budget=budget,
            suppliers=self._suppliers,
            max_products=_UNLIMITED,
            max_reactions=_UNLIMITED,
            max_compounds=max_compounds,
            max_iter=max_iter,
            shuffle=shuffle,
            debug=debug,
        )
        self._write_recipe(recipe, budget, stats, {'gen_max_compounds': max_compounds})
        return recipe


class RandomRecipeSelectionGenerator(_GeneratorBase):
    """Generate random recipes combining synthetic routes and compound selections."""

    def __init__(
        self,
        *,
        out_key: str,
        suppliers: 'list | None' = None,
        route_pool: 'RouteSet | None' = None,
        compounds=None,
        amount: float = 1.0,
        max_lead_time: 'float | None' = None,
        start_with: 'Recipe | None' = None,
    ) -> None:
        mrich.debug('RandomRecipeSelectionGenerator.__init__()')
        self._suppliers = suppliers
        self._max_lead_time = max_lead_time
        self._starting_recipe = start_with or Recipe()

        self._setup_io(
            out_key, '_rsgen.json', '_recipes_and_selections', skip_directory_creation=False
        )

        # inner generators build the pools and dump their own state files; they do
        # not create recipe directories (only this generator writes recipes)
        self._rgen = RandomRecipeGenerator(
            out_key=out_key,
            suppliers=suppliers,
            route_pool=route_pool,
            max_lead_time=max_lead_time,
            skip_directory_creation=True,
        )
        self._sgen = RandomSelectionGenerator(
            out_key=out_key,
            suppliers=suppliers,
            compounds=compounds,
            amount=amount,
            max_lead_time=max_lead_time,
            skip_directory_creation=True,
        )

        # combined pool: compound ingredients followed by routes
        self._pool = list(self._sgen.compound_pool) + list(self._rgen.route_pool)

        self._dump_data({})

    @property
    def compound_and_route_pool(self) -> list:
        """Combined pool of compound ingredients and routes"""
        return self._pool

    def generate(
        self,
        budget: float = 10000,
        currency: str = 'EUR',
        *,
        max_products: int = 1000,
        max_reactions: int = 1000,
        max_iter: int | None = None,
        shuffle: bool = True,
        debug: bool = False,
    ) -> 'Recipe':
        """Generate a random recipe of routes + compound selections (also written to disk)."""
        budget = Price(budget, currency)
        recipe, stats = _generate_recipe(
            self._starting_recipe,
            list(self._pool),
            budget=budget,
            suppliers=self._suppliers,
            max_products=max_products,
            max_reactions=max_reactions,
            max_compounds=_UNLIMITED,
            max_iter=max_iter,
            shuffle=shuffle,
            debug=debug,
        )
        self._write_recipe(
            recipe,
            budget,
            stats,
            {'gen_max_products': max_products, 'gen_max_reactions': max_reactions},
        )
        return recipe
