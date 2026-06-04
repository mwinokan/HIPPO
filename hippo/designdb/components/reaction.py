"""Reaction component."""

import mcol
import mrich
from designdb.models import CataloguePriceCompoundJunctionModel, CompoundModel, ReactionModel

from .compound import Compound

DEFAULT_REACTANT_AMOUNT = 1.0
DEFAULT_PRODUCT_YIELD = 1.0


class Reaction:
    """A :class:`.Reaction` wraps a :class:`.ReactionModel` and represents a
    synthetic step from reactant :class:`.Compound` objects to a product."""

    def __init__(self, instance: ReactionModel) -> None:
        self._instance = instance

    ### PROPERTIES

    @property
    def id(self) -> int:
        """Returns the :class:`.Reaction` ID"""
        return self._instance.pk

    @property
    def type(self) -> str:
        """Returns the reaction type string"""
        return self._instance.reaction_type

    @property
    def product_yield(self) -> float:
        """Returns the reaction product yield (fraction)"""
        return self._instance.reaction_product_yield

    @property
    def product(self) -> 'Compound':
        """Returns the product as a :class:`.Compound` component"""
        return Compound(self._instance.product_compound)

    @property
    def reactants(self) -> 'list[Compound]':
        """Returns the reactant :class:`.Compound` components"""
        return [Compound(r.compound) for r in self._instance.reactants.all()]

    @property
    def reactant_ids(self) -> list[int]:
        """Returns the reactant :class:`.CompoundModel` PKs"""
        return list(self._instance.reactants.values_list('compound_id', flat=True))

    @property
    def product_smiles(self) -> str:
        """Returns the product compound SMILES"""
        return self._instance.product_compound.compound_smiles

    @property
    def reaction_str(self) -> str:
        """Returns a human-readable reaction string"""
        s = ' + '.join(str(r) for r in self.reactants)
        return f'{s} -> {self.product}'

    @property
    def plain_repr(self) -> str:
        """Unformatted long string representation"""
        return f'{self}: {self.reaction_str} via {self.type}'

    ### METHODS

    def get_reactant_amount_pairs(self, compound_object: bool = True) -> list[tuple]:
        """Returns pairs of (compound, amount) for each reactant.

        :param compound_object: return :class:`.Compound` objects instead of IDs
        """
        pairs = [
            (cid, amount if amount is not None else DEFAULT_REACTANT_AMOUNT)
            for cid, amount in self._instance.reactants.values_list(
                'compound_id', 'reactant_amount'
            )
        ]
        if not pairs:
            return []
        if compound_object:
            return [
                (Compound(CompoundModel.objects.get(pk=cid)), amount)
                for cid, amount in pairs
            ]
        return pairs

    def check_reactant_availability(
        self,
        supplier: str | None = None,
        debug: bool = False,
    ) -> bool:
        """Check that every reactant either has a catalogue price or can be synthesised.

        :param supplier: restrict price check to this supplier
        :param debug: increase verbosity
        """
        for reactant_model in self._instance.reactants.all():
            compound = reactant_model.compound

            if debug:
                mrich.var('reactant', compound.pk)

            if supplier:
                has_quote = CataloguePriceCompoundJunctionModel.objects.filter(
                    compound=compound,
                    catalogue_price__supplier=supplier,
                ).exists()
            else:
                has_quote = CataloguePriceCompoundJunctionModel.objects.filter(
                    compound=compound,
                ).exists()

            has_reaction = ReactionModel.objects.filter(
                product_compound=compound
            ).exists()

            if debug:
                mrich.debug(f'{has_quote=}, {has_reaction=}')

            if not has_quote and not has_reaction:
                if debug:
                    mrich.warning(f'No quote or reaction for reactant pk={compound.pk}')
                return False

        return True

    def get_recipes(
        self,
        amount: float = 1,
        debug: bool = False,
        pick_cheapest: bool = False,
        permitted_reactions: 'ReactionSet | None' = None,
        supplier: str | None = None,
    ) -> 'Recipe | list[Recipe]':
        """Get a :class:`.Recipe` for this reaction.

        :param amount: amount in mg
        """
        # convenience bridge to the service layer (local import keeps the
        # component -> service dependency out of the module import graph)
        from designdb.services.recipe import RecipeService

        return RecipeService.from_reaction(
            self._instance,
            amount=amount,
            debug=debug,
            pick_cheapest=pick_cheapest,
            permitted_reactions=permitted_reactions,
            supplier=supplier,
        )

    ### DUNDERS

    def __str__(self) -> str:
        return f'R{self.id}'

    def __repr__(self) -> str:
        return (
            f'{mcol.bold}{mcol.underline}{self.plain_repr}'
            f'{mcol.unbold}{mcol.ununderline}'
        )

    def __rich__(self) -> str:
        return f'[bold underline]{self.plain_repr}'

    def __eq__(self, other: 'int | Reaction | ReactionModel') -> bool:
        match other:
            case int():
                return self.id == other
            case Reaction():
                return self._instance == other._instance
            case ReactionModel():
                return self._instance == other
            case _:
                raise NotImplementedError

    def __hash__(self) -> int:
        return self.id
