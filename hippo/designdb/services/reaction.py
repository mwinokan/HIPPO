import logging

import mrich

# from mypackage.services.compound import CompoundService
# from rdkit.Chem import inchi
from designdb.models import CompoundModel, ReactantModel, ReactionModel

logger = logging.getLogger(__name__)


class ReactionService:
    @staticmethod
    def reactant_compound_ids() -> set[int]:
        """Return compound IDs that are reactants of at least one reaction and not
        a product of any (i.e. leaf reactants / purchasable building blocks).

        Mirrors the legacy ``CompoundSet.reactants`` (reactant compounds minus
        compounds that are a reaction product).
        """
        reactant_ids = set(
            ReactantModel.objects.values_list('compound_id', flat=True).distinct()
        )
        product_ids = set(
            ReactionModel.objects.values_list(
                'product_compound_id', flat=True
            ).distinct()
        )
        return reactant_ids - product_ids

    @classmethod
    def create_from_lists(
        cls,
        *,
        reaction_types: list[str],
        product_ids: list[int],
        reactant_id_lists: list[set[int]],
    ) -> list[int]:
        # insert reaction

        # insert reactant

        reaction_ids = []
        non_duplicates = {}

        # not entirely sure how the original query was meant to work
        qs = ReactantModel.objects.filter(compound__pk__in=product_ids)
        existing = {}
        for r in qs:
            reaction_type = r.reaction.reaction_type
            reaction_product = r.reaction.product_compound.pk
            reaction_id = r.reaction.pk
            reactant_compound = r.compound.pk

            key = (reaction_type, reaction_product)

            if key not in existing:
                existing[key] = {}

            if reaction_id not in existing[key]:
                existing[key][reaction_id] = set()

            existing[key][reaction_id].add(reactant_compound)

        existing_count = 0

        # why is strict false??
        for reaction_type, product_id, reactant_ids in zip(
            reaction_types, product_ids, reactant_id_lists, strict=False
        ):
            key = (reaction_type, product_id)

            possible_matches = {k: v for k, v in existing.items() if k == key}

            assert len(possible_matches) < 2

            if possible_matches:
                possible_matches = list(possible_matches.values())[0]

            if any(reactant_ids == v for v in possible_matches.values()):
                existing_count += 1
                continue

            non_duplicates[key] = reactant_ids

        if existing_count:
            mrich.warning('Skipped', existing_count, 'existing reactions')

        if not non_duplicates:
            mrich.warning('All reactions are duplicates')
            return None

        for reaction_type, product_id in non_duplicates.keys():
            compound = CompoundModel.objects.get(pk=product_id)
            # if I understand the original procedure correctly, it
            # should have already weeded out the duplicates
            reaction, _ = ReactionModel.objects.get_or_create(
                reaction_type=reaction_type,
                product_compound=compound,
                reaction_product_yield=1.0,
            )
            reaction_ids.append(reaction.pk)

        payload = []
        for reaction_id, ((_reaction_type, _product_id), reactant_ids) in zip(
            reaction_ids, non_duplicates.items(), strict=False
        ):
            for reactant_id in reactant_ids:
                payload.append((reaction_id, reactant_id))

        for reaction_id, reactant_id in payload:
            reaction = ReactionModel.objects.get(pk=reaction_id)
            compound = CompoundModel.objects.get(pk=reactant_id)
            reaction, _ = ReactantModel.objects.get_or_create(
                reaction=reaction,
                compound=compound,
                reactant_amount=1.0,
            )

        # delete orphaned reactions, srsly??
        ReactionModel.objects.filter(
            pk__in=ReactantModel.objects.filter(
                compound__isnull=True,
            ).values('reaction'),
        ).delete()

        return reaction_ids
