# from mypackage.services.compound import CompoundService

# from rdkit.Chem import inchi
import logging
from collections import Counter

import mrich
from designdb.models import ComponentModel, RouteModel
from designdb.recipe import Recipe

logger = logging.getLogger(__name__)


class RouteService:
    @classmethod
    def create_from_recipe(
        cls,
        *,
        recipe: Recipe,
    ) -> tuple[RouteModel, bool]:

        route, created = RouteModel.objects.get_or_create(
            product_compound=recipe.product.compound
        )

        # are you joking?? reactants and intermediates are all of the
        # sudden components

        # component_type encoding (see Route.get_route): 1=reaction, 2=reactant,
        # 3=intermediate
        components = []

        # reactions
        components.extend(
            [
                ComponentModel(route=route, component_type=1, component_ref=ref.pk)
                for ref in recipe.reactions
            ],
        )

        # reactants
        components.extend(
            [
                ComponentModel(
                    route=route,
                    component_type=2,
                    component_ref=ref,
                    component_amount=amount,
                )
                for ref, amount in recipe.reactants.id_amount_pairs
            ],
        )

        # intermediates
        components.extend(
            [
                ComponentModel(
                    route=route,
                    component_type=3,
                    component_ref=ref,
                    component_amount=amount,
                )
                for ref, amount in recipe.intermediates.id_amount_pairs
            ],
        )

        ComponentModel.objects.bulk_create(components, ignore_conflicts=True)

        return route, created

    @classmethod
    def prune_duplicate_routes(cls) -> int:
        """Delete duplicate routes, keeping the lowest-id copy of each.

        A duplicate is any pair of routes with the same product compound and
        identical sets of (component_ref, component_type) pairs.

        Returns the number of routes deleted.
        """
        rows = ComponentModel.objects.values_list(
            'route_id', 'route__product_compound_id', 'component_ref', 'component_type'
        )

        route_fingerprints: dict[int, tuple] = {}
        for route_id, product_id, comp_ref, comp_type in rows:
            if route_id not in route_fingerprints:
                route_fingerprints[route_id] = (product_id, set())
            route_fingerprints[route_id][1].add((comp_ref, comp_type))

        # freeze the sets so they're hashable
        frozen = {rid: (fp[0], frozenset(fp[1])) for rid, fp in route_fingerprints.items()}

        mrich.var('#routes', len(frozen))

        counter = Counter(frozen.values())
        duplicates = {fp: count for fp, count in counter.items() if count > 1}
        mrich.var('products with duplicate routes', len(duplicates))

        if not duplicates:
            mrich.success('No duplicate routes found')
            return 0

        to_delete: set[int] = set()
        for fp in duplicates:
            matched = sorted(rid for rid, v in frozen.items() if v == fp)
            mrich.print('compound', fp[0], 'has', len(matched), 'duplicate routes')
            to_delete.update(matched[1:])

        ComponentModel.objects.filter(route_id__in=to_delete).delete()
        deleted, _ = RouteModel.objects.filter(pk__in=to_delete).delete()
        mrich.success('Deleted', len(to_delete), 'duplicate routes')
        return len(to_delete)

    # @property
    # def id_amount_pairs(self) -> list[tuple]:
    #     """Get a list of compound ID and amount pairs"""
    #     return [
    #         (id, amount) for id, amount in self.df[['compound_id', 'amount']].values
    #     ]
