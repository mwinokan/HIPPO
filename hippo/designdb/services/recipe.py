"""Recipe construction and DB traversal.

Builds :class:`.Recipe` objects from reactions/compounds/reactants. The
:class:`.Recipe` aggregate itself is lean; its ``from_*`` classmethods are
deprecated shims that delegate here.
"""

from itertools import product
from typing import TYPE_CHECKING

import mrich
from designdb.components.compound import Compound
from designdb.components.reaction import DEFAULT_PRODUCT_YIELD, Reaction
from designdb.models import (
    CompoundModel,
    InspirationModel,
    PoseModel,
    ReactionModel,
    RouteModel,
)
from designdb.sets.compound import CompoundSet
from designdb.sets.ingredient import IngredientSet
from designdb.sets.pose import PoseSet
from designdb.sets.reaction import ReactionSet

if TYPE_CHECKING:
    from pathlib import Path

    from designdb.recipe import Recipe
    from designdb.sets.route import RouteSet
    from pandas import DataFrame


class RecipeService:
    """Construction and traversal logic for :class:`.Recipe` objects."""

    ### FACTORIES

    @staticmethod
    def from_reaction(
        reaction,
        amount=1,
        *,
        debug: bool = False,
        pick_cheapest: bool = True,
        permitted_reactions: 'ReactionSet | None' = None,
        quoted_only: bool = False,
        supplier: None | str = None,
        unavailable_reaction: str = 'error',
        reaction_checking_cache: dict[int, bool] | None = None,
        reaction_reactant_cache: dict[int, bool] | None = None,
        inner: bool = False,
        get_ingredient_quotes: bool = True,
    ) -> 'Recipe | list[Recipe] | None':
        """Create a :class:`.Recipe` from a :class:`.ReactionModel` and its upstream
        dependencies.

        :param reaction: :class:`.ReactionModel` to create the recipe from
        :param amount: amount in ``mg`` (Default value = 1)
        :param debug: increase verbosity (Default value = False)
        :param pick_cheapest: return only the cheapest solution (Default value = True)
        :param permitted_reactions: only consider reactions in this set
        :param quoted_only: only allow reactants with quotes (Default value = False)
        :param supplier: restrict quotes to this supplier (Default value = None)
        :param unavailable_reaction: behaviour when a reaction has unavailable
            reactants (Default value = 'error')
        :param inner: indicates a recursive call (Default value = False)
        :param get_ingredient_quotes: get quotes for product ingredients
        """

        from designdb.recipe import Recipe

        assert isinstance(reaction, ReactionModel)
        reaction_c = Reaction(reaction)

        if debug:
            mrich.debug(
                f'RecipeService.from_reaction(R{reaction.id}, '
                f'{amount=}, {pick_cheapest=})'
            )
            mrich.debug(f'{reaction_c.product.id=}')
            mrich.debug(f'{reaction_c.reactant_ids=}')

        if permitted_reactions:
            assert reaction in permitted_reactions

        recipe = Recipe(
            products=IngredientSet(
                [
                    reaction_c.product.as_ingredient(
                        amount=amount, get_quote=get_ingredient_quotes
                    )
                ],
            ),
            reactants=IngredientSet([], supplier=supplier),
            intermediates=IngredientSet([]),
            reactions=ReactionSet([reaction.id], sort=False),
        )

        recipes = [recipe]

        if quoted_only or supplier:
            if debug:
                mrich.debug(f'Checking reactant_availability: {reaction=}')
            if reaction_checking_cache and reaction.id in reaction_checking_cache:
                ok = reaction_checking_cache[reaction.id]
            else:
                ok = reaction_c.check_reactant_availability(supplier=supplier)
                if reaction_checking_cache is not None:
                    reaction_checking_cache[reaction.id] = ok
            if not ok:
                if unavailable_reaction == 'error':
                    mrich.error(f'Reactants not available for {reaction=}')
                return None if pick_cheapest else []

        def get_reactant_amount_pairs(
            reaction_model: 'ReactionModel',
        ) -> list[tuple[int, float]]:
            """Get pairs of reactant ID and float amounts"""
            if reaction_reactant_cache and reaction_model.id in reaction_reactant_cache:
                return reaction_reactant_cache[reaction_model.id]
            pairs = Reaction(reaction_model).get_reactant_amount_pairs(
                compound_object=False
            )
            if reaction_reactant_cache is not None:
                reaction_reactant_cache[reaction_model.id] = pairs
            return pairs

        if debug:
            mrich.debug(f'get_reactant_amount_pairs({reaction.id})')
        pairs = get_reactant_amount_pairs(reaction)

        for reactant_id, reactant_amount in pairs:
            reactant = Compound(CompoundModel.objects.get(pk=reactant_id))

            if debug:
                mrich.debug(f'{reactant.id=}, {reactant_amount=}')

            # scale amount
            reactant_amount *= amount
            reactant_amount /= reaction_c.product_yield or DEFAULT_PRODUCT_YIELD

            inner_reactions = reactant.get_reactions(
                none='quiet', permitted_reactions=permitted_reactions
            )

            if len(inner_reactions):
                if debug:
                    if len(inner_reactions) == 1:
                        mrich.debug('Reactant has ONE inner reaction')
                    else:
                        mrich.warning(f'{reactant=} has MULTIPLE inner reactions')

                inner_recipes = []
                for inner_reaction in inner_reactions:
                    reaction_recipes = RecipeService.from_reaction(
                        reaction=inner_reaction,
                        amount=reactant_amount,
                        debug=debug,
                        pick_cheapest=False,
                        quoted_only=quoted_only,
                        supplier=supplier,
                        unavailable_reaction=unavailable_reaction,
                        reaction_checking_cache=reaction_checking_cache,
                        reaction_reactant_cache=reaction_reactant_cache,
                        inner=True,
                    )
                    inner_recipes += reaction_recipes

                new_recipes = []
                for recipe in recipes:
                    for inner_recipe in inner_recipes:
                        combined_recipe = recipe.copy()

                        combined_recipe.reactants += inner_recipe.reactants
                        combined_recipe.intermediates += inner_recipe.intermediates
                        combined_recipe.reactions += inner_recipe.reactions
                        combined_recipe.intermediates.add(
                            reactant.as_ingredient(reactant_amount, supplier=supplier)
                        )

                        new_recipes.append(combined_recipe)

                recipes = new_recipes

            else:
                ingredient = reactant.as_ingredient(reactant_amount, supplier=supplier)
                for recipe in recipes:
                    recipe.reactants.add(ingredient)

        # reverse ReactionSet's (outermost call only)
        if not inner:
            for recipe in recipes:
                recipe.reactions.reverse()

        if pick_cheapest:
            if debug:
                mrich.debug('Picking cheapest')
            priced = [r for r in recipes if r.get_price(supplier=supplier)]
            if not priced:
                mrich.error("0 recipes with prices, can't choose cheapest")
                return recipes
            sorted_recipes = sorted(
                priced, key=lambda r: r.get_price(supplier=supplier)
            )
            if debug:
                for recipe in recipes:
                    mrich.debug(f'{recipe}, {recipe.price}')
            return sorted_recipes[0]

        return recipes

    @staticmethod
    def from_reactions(
        reactions: 'ReactionSet',
        amount: float = 1,
        pick_cheapest: bool = True,
        permitted_reactions: 'ReactionSet | None' = None,
        final_products_only: bool = True,
        return_products: bool = False,
        supplier: str | None = None,
        use_routes: bool = False,
        debug: bool = False,
        **kwargs,
    ) -> 'Recipe | list[Recipe] | CompoundSet':
        """Create a :class:`.Recipe` from a :class:`.ReactionSet` and its upstream
        dependencies.

        :param reactions: reactions to create the recipe from
        :param amount: amount in ``mg`` (Default value = 1)
        :param pick_cheapest: choose the cheapest solution (Default value = True)
        :param permitted_reactions: only consider reactions in this set
        :param final_products_only: don't make routes to intermediates
            (Default value = True)
        :param return_products: return the :class:`.CompoundSet` of products instead
        """

        assert isinstance(reactions, ReactionSet)

        if debug:
            mrich.debug('RecipeService.from_reactions()')
            mrich.var('reactions', reactions)
            mrich.var('amount', amount)
            mrich.var('final_products_only', final_products_only)
            mrich.var('permitted_reactions', permitted_reactions)

        # all products synthesisable from these reactions
        products = reactions.products

        if debug:
            mrich.var('products', products)

        if final_products_only:
            # keep only compounds that are never used as a reactant (i.e. leaves)
            from designdb.models import ReactantModel

            products = CompoundSet(
                CompoundModel.objects.filter(pk__in=list(products.ids)).exclude(
                    pk__in=ReactantModel.objects.values('compound'),
                )
            )
            if debug:
                mrich.var('final products', products)

            if return_products:
                return products

        return RecipeService.from_compounds(
            compounds=products,
            amount=amount,
            permitted_reactions=reactions,
            pick_cheapest=pick_cheapest,
            supplier=supplier,
            use_routes=use_routes,
            debug=debug,
            **kwargs,
        )

    @staticmethod
    def from_compounds(
        compounds: 'CompoundSet',
        amount: float = 1,
        debug: bool = False,
        pick_cheapest: bool = True,
        permitted_reactions: 'ReactionSet | None' = None,
        quoted_only: bool = False,
        supplier: None | str = None,
        solve_combinations: bool = True,
        pick_first: bool = False,
        warn_multiple_solutions: bool = True,
        pick_cheapest_inner_routes: bool = False,
        unavailable_reaction: str = 'error',
        reaction_checking_cache: dict[int, bool] | None = None,
        reaction_reactant_cache: dict[int, bool] | None = None,
        use_routes: bool = False,
        **kwargs,
    ):
        """Create recipe(s) to synthesise the products in a :class:`.CompoundSet`.

        :param compounds: set of compounds to find routes for
        :param solve_combinations: combinatorially combine the individual solutions
            (Default value = True)
        :param pick_first: return the first solution without comparison
        :param warn_multiple_solutions: warn if a compound has multiple routes
        :param pick_cheapest_inner_routes: for each compound choose the cheapest route
        :param use_routes: use stored :class:`.RouteModel` rows instead of solving
            reactions on the fly
        """

        from designdb.recipe import Route

        assert isinstance(compounds, CompoundSet)

        n_comps = len(compounds)
        assert n_comps

        if not hasattr(amount, '__iter__'):
            amount = [amount] * n_comps

        if use_routes and supplier:
            raise NotImplementedError(
                'use_routes combined with a supplier filter is not supported'
            )

        options = []
        ok = 0
        mrich.var('#compounds', n_comps)

        for comp, a in mrich.track(
            zip(compounds, amount, strict=False),
            prefix='Solving individual compound recipes...',
            total=n_comps,
        ):
            comp_options = []

            if use_routes:
                route_ids = list(
                    RouteModel.objects.filter(product_compound__id=comp.id).values_list(
                        'id', flat=True
                    )
                )
                if not route_ids:
                    mrich.error('No routes to', comp)
                    continue
                comp_options = [Route.get_route(id=route_id) for route_id in route_ids]

            else:
                for reaction in Compound(comp).reactions:
                    if permitted_reactions and reaction not in permitted_reactions:
                        continue

                    sol = RecipeService.from_reaction(
                        reaction=reaction,
                        amount=a,
                        pick_cheapest=pick_cheapest_inner_routes,
                        debug=debug,
                        permitted_reactions=permitted_reactions,
                        quoted_only=quoted_only,
                        supplier=supplier,
                        unavailable_reaction=unavailable_reaction,
                        reaction_checking_cache=reaction_checking_cache,
                        reaction_reactant_cache=reaction_reactant_cache,
                        **kwargs,
                    )

                    if pick_cheapest_inner_routes:
                        if sol:
                            comp_options.append(sol)
                    else:
                        assert isinstance(sol, list)
                        comp_options += sol

                if not comp_options:
                    mrich.error(
                        f'No solutions for compound={comp} '
                        f'({Compound(comp).reactions.ids=})'
                    )
                    continue

            if pick_cheapest and len(comp_options) > 1:
                if warn_multiple_solutions:
                    mrich.warning(
                        'Multiple solutions for', comp, '(', len(comp_options), ')'
                    )
                if debug:
                    mrich.debug('Picking cheapest...')
                priced = [r for r in comp_options if r.price]
                comp_options = sorted(priced, key=lambda r: r.price)[:1]

            if warn_multiple_solutions and len(comp_options) > 1:
                mrich.warning(f'Multiple solutions for compound={comp}')
                if debug:
                    mrich.debug(f'{comp_options=}')
            else:
                if n_comps <= 200:
                    mrich.success(f'Found solution for compound={comp}')
                ok += 1
                mrich.set_progress_field('ok', ok)
                mrich.set_progress_field('n', n_comps)

            options.append(comp_options)

        assert all(options)

        mrich.print('Solving recipe combinations...')
        combinations = list(product(*options))

        if not solve_combinations:
            return combinations

        solutions = []

        if n_comps > 1:
            generator = mrich.track(
                combinations, prefix='Combining recipes...', total=len(combinations)
            )
        else:
            generator = combinations

        ok = 0
        for combo in generator:
            if debug:
                mrich.debug(f'Combination of {len(combo)} recipes')

            if not combo:
                continue

            solution = combo[0]
            for i, recipe in enumerate(combo[1:]):
                if debug:
                    mrich.debug(i + 1)
                solution += recipe

            solutions.append(solution)
            ok += 1
            mrich.set_progress_field('ok', ok)
            mrich.set_progress_field('n', len(combinations))

        if not solutions:
            mrich.error('No solutions')
            return None

        if pick_first:
            return solutions[0]

        if pick_cheapest:
            mrich.debug('Calculating prices...')
            priced = [r for r in solutions if r.price]
            mrich.print('Picking cheapest from', len(priced), 'options')
            if not priced:
                mrich.error("0 recipes with prices, can't choose cheapest")
                # fall back to the first (unpriced) solution so the return type
                # stays a single Recipe, consistent with pick_first / the priced path
                return solutions[0]
            return sorted(priced, key=lambda r: r.price)[0]

        return solutions

    @staticmethod
    def from_reactants(
        reactants: 'CompoundSet | IngredientSet',
        amount: float = 1,
        debug: bool = False,
        return_products: bool = False,
        supplier: str | None = None,
        pick_cheapest: bool = False,
        use_routes: bool = False,
        **kwargs,
    ) -> 'list[Recipe] | Recipe | CompoundSet':
        """Find the maximal recipe reachable from a given set of reactants.

        :param reactants: :class:`.CompoundSet` or :class:`.IngredientSet` of
            reactants (ingredient amounts are ignored)
        :param amount: amount of each product needed (Default value = 1)
        :param return_products: return the products instead of the recipe
        """

        if isinstance(reactants, IngredientSet):
            reactant_ids = reactants.compound_ids
        else:
            reactant_ids = reactants.ids

        all_reactants = set(reactant_ids)
        possible_reactions: set[int] = set()

        # recursively expand the set of reachable reactions/products
        for _ in range(300):
            reaction_ids = RecipeService._possible_reaction_ids(all_reactants)

            if not reaction_ids:
                break

            if debug:
                mrich.debug(f'Adding {len(reaction_ids)} reactions')

            possible_reactions |= set(reaction_ids)

            product_ids = list(
                ReactionModel.objects.filter(pk__in=reaction_ids).values_list(
                    'product_compound_id', flat=True
                )
            )

            n_prev = len(all_reactants)
            all_reactants |= set(product_ids)

            if n_prev == len(all_reactants):
                break
        else:
            raise NotImplementedError('Maximum recursion depth exceeded')

        if debug:
            mrich.var('all possible reactions', possible_reactions)

        rset = ReactionSet(list(possible_reactions), sort=False)

        return RecipeService.from_reactions(
            rset,
            amount=amount,
            permitted_reactions=rset,
            debug=debug,
            return_products=return_products,
            supplier=supplier,
            use_routes=use_routes,
            pick_cheapest=pick_cheapest,
            **kwargs,
        )

    ### TRAVERSAL

    @staticmethod
    def get_routes(recipe: 'Recipe', return_ids: bool = False) -> 'RouteSet':
        """Get stored routes to the products of ``recipe`` restricted to its
        reactions."""
        return recipe.products.compounds.get_routes(
            permitted_reactions=recipe.reactions, return_ids=return_ids
        )

    ### EXPORTERS

    @staticmethod
    def write_CAR_csv(
        recipe: 'Recipe', file: 'str | Path', return_df: bool = False
    ) -> 'DataFrame | None':
        """Write CSV(s) for use with CAR.

        Requires a populated ``route`` table (see :meth:`.RecipeService.get_routes`).
        One row per route; reactions are flattened into ``reactant-N-i`` /
        ``reaction-product-smiles-i`` / ``reaction-name-i`` columns.

        :param recipe: the :class:`.Recipe` to export
        :param file: output path (per-step files are also written alongside)
        :param return_df: return the assembled DataFrame
        """

        from pathlib import Path

        from pandas import DataFrame

        file = str(Path(file).resolve())
        rows = []

        for sub_recipe in RecipeService.get_routes(recipe):
            product = sub_recipe.product

            row = {
                'target-names': str(product.compound),
                'no-steps': 0,
                'concentration-required-mM': None,
                'amount-required-uL': None,
                'batch-tag': None,
            }

            for i, reaction_model in enumerate(sub_recipe.reactions):
                i = i + 1
                reaction = Reaction(reaction_model)
                row['no-steps'] += 1

                reactants = reaction.reactants
                match len(reactants):
                    case 1:
                        row[f'reactant-1-{i}'] = reactants[0].smiles
                        row[f'reactant-2-{i}'] = None
                    case 2:
                        row[f'reactant-1-{i}'] = reactants[0].smiles
                        row[f'reactant-2-{i}'] = reactants[1].smiles
                    case _:
                        for j, reactant in enumerate(reactants):
                            row[f'reactant-{j + 1}-{i}'] = reactant.smiles

                row[f'reaction-product-smiles-{i}'] = reaction.product_smiles
                row[f'reaction-name-{i}'] = reaction.type
                row[f'reaction-recipe-{i}'] = None
                row[f'reaction-groupby-column-{i}'] = None

            rows.append(row)

        df = DataFrame(rows)

        if len(df[df.duplicated()]):
            mrich.warning('Removing duplicates from CAR DataFrame')
            df = df.drop_duplicates()

        df = df.convert_dtypes()

        for n_steps in set(df['no-steps']):
            subset = df[df['no-steps'] == n_steps]
            this_file = file.replace('.csv', f'_{n_steps}steps.csv')
            mrich.writing(this_file)
            subset.to_csv(this_file, index=False)

        mrich.writing(file)
        df.to_csv(file, index=False)

        if return_df:
            return df
        return None

    @staticmethod
    def write_reactant_csv(recipe: 'Recipe', file, reaction_type_counts=True, **kwargs):
        """Detailed reactant-purchasing CSV. Not implemented."""
        raise NotImplementedError('write_reactant_csv is not implemented')

    @staticmethod
    def write_product_csv(
        recipe: 'Recipe', file, return_df: bool = False
    ) -> 'DataFrame | None':
        """Detailed CSV output including product information for selection/synthesis.

        One row per product compound: identifiers, required amount, associated poses,
        tags, upstream route/reaction/reactant dependencies, scaffold series, and
        inspiration pose names.

        :param recipe: the :class:`.Recipe` whose products to report
        :param file: output CSV path
        :param return_df: also return the assembled ``DataFrame``
        """

        from pandas import DataFrame

        routes = RecipeService.get_routes(recipe)

        product_ids = list(recipe.products.compound_ids)

        # compound_id -> set of associated pose IDs
        pose_map: dict[int, set] = {}
        for comp_id, pose_id in PoseModel.objects.filter(
            compound_id__in=product_ids
        ).values_list('compound_id', 'id'):
            pose_map.setdefault(comp_id, set()).add(pose_id)

        # compound_id -> set of inspiration (original) pose IDs, scoped to the
        # product compounds and their scaffolds (the inspiration fallback needs both)
        scaffold_ids = set()
        for prod in recipe.products:
            # Ingredient.__getattr__ delegates to the CompoundModel (ORM), so wrap
            # in the Compound component to reach component-level properties
            if scaffolds := Compound(prod.compound).scaffolds:
                scaffold_ids.update(scaffolds.ids)
        needed_ids = set(product_ids) | scaffold_ids

        inspiration_map: dict[int, set] = {}
        for comp_id, original_pose_id in InspirationModel.objects.filter(
            derivative_pose__compound_id__in=needed_ids
        ).values_list('derivative_pose__compound_id', 'original_pose_id'):
            inspiration_map.setdefault(comp_id, set()).add(original_pose_id)

        data = []

        for prod in mrich.track(
            recipe.products, prefix='Constructing product DataFrame'
        ):
            # wrap in the Compound component for component-level properties
            # (Ingredient.__getattr__ delegates to the CompoundModel ORM instead)
            comp = Compound(prod.compound)

            d = dict(
                hippo_id=prod.compound_id,
                smiles=comp.smiles,
                inchikey=comp.inchikey,
                required_amount_mg=prod.amount,
            )

            upstream_routes = []
            upstream_reaction_ids = []

            for route in routes:
                if route.product_compound.id == prod.compound_id:
                    upstream_routes.append(route)
                    upstream_reaction_ids += route.reactions.ids

            if not upstream_routes:
                mrich.error('No upstream routes for', prod)
                continue

            if not upstream_reaction_ids:
                mrich.error('No upstream reactions for', prod)
                continue

            upstream_reactions = ReactionSet(list(set(upstream_reaction_ids)))

            # scaffold series: the product's scaffolds, or itself if it is one
            if scaffolds := comp.scaffolds:
                scaffold_series, is_scaffold = scaffolds.ids, False
            else:
                scaffold_series, is_scaffold = [prod.compound_id], True

            poses = pose_map.get(prod.compound_id, set())

            d['num_poses'] = len(poses)
            d['poses'] = poses
            d['tags'] = comp.tags
            d['num_routes'] = len(upstream_routes)
            d['num_reaction_steps'] = {len(r.reactions) for r in upstream_routes}
            d['reaction_dependencies'] = upstream_reactions.ids
            d['reactant_dependencies'] = set(
                sum((route.reactants.ids for route in upstream_routes), [])
            )
            d['route_ids'] = [route.id for route in upstream_routes]
            d['chemistry_types'] = ', '.join(t for t in upstream_reactions.types if t)
            d['is_scaffold'] = is_scaffold
            d['scaffold_series'] = scaffold_series

            # inspiration pose IDs, with fallback to the scaffold / metadata
            inspirations = inspiration_map.get(prod.compound_id, None)

            if not inspirations and not is_scaffold:
                scaffold = Compound(comp.scaffolds[0])
                inspirations = inspiration_map.get(scaffold.id, None)

                scaffold_meta = scaffold.metadata or {}
                if not inspirations and 'inspiration_pose_ids' in scaffold_meta:
                    inspirations = scaffold_meta['inspiration_pose_ids']

            product_meta = comp.metadata or {}
            if (
                not inspirations
                and is_scaffold
                and 'inspiration_pose_ids' in product_meta
            ):
                inspirations = product_meta['inspiration_pose_ids']

            if inspirations:
                inspiration_poses = PoseSet(
                    PoseModel.objects.filter(pk__in=list(inspirations))
                )
                d['inspirations'] = ', '.join(inspiration_poses.names)
            else:
                d['inspirations'] = ''

            data.append(d)

        df = DataFrame(data)
        mrich.writing(file)
        df.to_csv(file, index=False)

        if return_df:
            return df
        return None

    @staticmethod
    def to_syndirella(recipe: 'Recipe', out_key, poses, *, separate: bool = False):
        """Generate Syndirella elaboration inputs from this recipe. Not implemented."""
        raise NotImplementedError('RecipeService.to_syndirella is not implemented')

    @staticmethod
    def register_missing_routes(
        recipe: 'Recipe', missing_only: bool = True, supplier: str = 'Enamine'
    ) -> None:
        """Calculate and register missing routes to the products of ``recipe``.
        Not implemented."""
        raise NotImplementedError('register_missing_routes is not implemented')

    ### HELPERS

    @staticmethod
    def _possible_reaction_ids(compound_ids: set[int]) -> list[int]:
        """Return reaction IDs whose every reactant is in ``compound_ids``."""
        from designdb.models import ReactantModel

        compound_ids = set(compound_ids)

        # reactions that use at least one of these compounds as a reactant
        candidate_ids = (
            ReactantModel.objects.filter(compound_id__in=compound_ids)
            .values_list('reaction_id', flat=True)
            .distinct()
        )

        # of those, keep reactions whose reactants are all available
        rows = ReactantModel.objects.filter(
            reaction_id__in=list(candidate_ids)
        ).values_list('reaction_id', 'compound_id')

        reaction_reactants: dict[int, set[int]] = {}
        for reaction_id, compound_id in rows:
            reaction_reactants.setdefault(reaction_id, set()).add(compound_id)

        return [
            reaction_id
            for reaction_id, reactants in reaction_reactants.items()
            if reactants <= compound_ids
        ]
