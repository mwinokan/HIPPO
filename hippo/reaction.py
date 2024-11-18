import mcol

from .compound import Compound
from .recipe import Recipe

import mrich


class Reaction:
    """
    A :class:`.Reaction` is a simplified representation of a synthetic pathway to create a product :class:`.Compound`. Reactants (also :class:`.Compound` objects) as well as a reaction type are required.

    .. attention::

            :class:`.Reaction` objects should not be created directly. Instead use :meth:`.HIPPO.register_reaction` or :meth:`.HIPPO.reactions`

    """

    _table = "reaction"

    def __init__(
        self,
        db: "Database",
        id: int,
        type: str,
        product: int,
        product_yield: float,
    ) -> None:

        self._db = db
        self._id = id
        self._type = type
        self._product_id = product
        self._product = None
        self._product_yield = product_yield

    ### PROPERTIES

    @property
    def id(self) -> int:
        """Returns the :class:`.Reaction` ID"""
        return self._id

    @property
    def type(self) -> str:
        """Returns the :class:`.Reaction` tyoe"""
        return self._type

    @property
    def product(self) -> "Compound":
        """Returns the reaction's product :class:`.Compound`"""
        if self._product is None:
            self._product = self.db.get_compound(id=self.product_id)
        return self._product

    @property
    def product_yield(self) -> float:
        """Returns the reaction's product yield (fraction)"""
        return self._product_yield

    @property
    def db(self) -> "Database":
        """Returns a pointer to the parent database"""
        return self._db

    @property
    def reactants(self) -> "CompoundSet":
        """Returns a :class:`.CompoundSet` of the reactants"""
        from .cset import CompoundSet

        return CompoundSet(self.db, indices=self.reactant_ids)

    @property
    def reaction_str(self) -> str:
        """Returns a string representing the reaction"""
        s = " + ".join([str(r) for r in self.reactants])
        s = f"{s} -> {str(self.product)}"
        return s

    @property
    def reactant_ids(self) -> set[int]:
        """Returns a set of reactant ID's"""
        return set(v for v in self.get_reactant_ids())

    @property
    def reactant_str_ids(self) -> str:
        """Return an SQL formatted tuple string of the reactant :class:`.Compound` IDs"""
        return str(tuple(self.reactant_ids)).replace(",)", ")")

    @property
    def product_id(self) -> int:
        """Returns the product :class:`.Compound` ID"""
        return self._product_id

    @property
    def product_smiles(self) -> str:
        """Product :class:`.Compound` SMILES string"""
        return self.product.smiles

    @property
    def reactant_smiles(self) -> list[str]:
        """List of reactant :class:`.Compound` SMILES strings"""
        return [r.smiles for r in self.reactants]

    @property
    def product_mol(self):
        """Product :class:`.Compound` ``rdkit.Chem.Mol`` object"""
        return self.product.mol

    @property
    def reactant_mols(self):
        """List of reactant :class:`.Compound` ``rdkit.Chem.Mol`` object"""
        return [r.mol for r in self.reactants]

    @property
    def price_estimate(self) -> float:
        """Estimate the price of this :class:`.Reaction`"""
        return self.db.get_reaction_price_estimate(reaction=self)

    @property
    def plain_repr(self) -> str:
        """Unformatted long string representation"""
        return f"{self}: {self.reaction_str} via {self.type}"

    ### METHODS

    def get_reactant_amount_pairs(self, compound_object: bool = True) -> list[tuple]:
        """Returns pairs of reactants and their amounts

        :param compound_object: return :class:`.Compound` object instead of ID, (Default value = True)
        :returns: list of tuples containing :class:`.Compound` ID/object and amount in mg
        """

        compound_ids = self.db.select_where(
            query="reactant_compound, reactant_amount",
            table="reactant",
            key="reaction",
            value=self.id,
            multiple=True,
        )

        if compound_ids:
            if compound_object:
                return [
                    # (self.db.get_compound(id=id), amount/self.product_yield) for id, amount in compound_ids
                    (self.db.get_compound(id=id), amount)
                    for id, amount in compound_ids
                ]
            else:
                return compound_ids
        else:
            return []

    def get_reactant_ids(self) -> list[int]:
        """Returns list of reactants :class:`.Compound` IDs

        :returns: list of :class:`.Compound` IDs
        """

        compound_ids = self.db.select_where(
            query="reactant_compound",
            table="reactant",
            key="reaction",
            value=self.id,
            multiple=True,
        )

        if compound_ids:
            return [id for id, in compound_ids]
        else:
            return []

    def get_recipes(
        self,
        amount: float = 1,  # in mg
        debug: bool = False,
        pick_cheapest: bool = False,
        permitted_reactions: "None | ReactionSet" = None,
        supplier: str | None = None,
    ) -> "Recipe | list[Recipe]":
        """Get a :class:`.Recipe` describing how to make the product

        :param amount: Amount in ``mg``, defaults to ``1``
        :param debug: Increase verbosity, (Default value = False)
        :param pick_cheapest: pick the cheapest :class:`.Recipe`, (Default value = False)
        :param permitted_reactions: Limit the reactions to consider to members of this set, (Default value = None)
        :param supplier: Limit to reactants from this supplier (Default value = None)
        :returns: :class:`.Recipe` object or list thereof

        """

        from .recipe import Recipe

        return Recipe.from_reaction(
            self,
            amount=amount,
            debug=debug,
            pick_cheapest=pick_cheapest,
            permitted_reactions=permitted_reactions,
            supplier=supplier,
        )

    def summary(
        self,
        draw: bool = True,
    ) -> None:
        """Print a summary of this reaction's information

        :param draw: draw the reaction compounds (Default value = True)

        """

        print(f"id={self.id}")
        print(f"type={self.type}")
        print(f"product={self.product}")
        print(f"product_yield={self.product_yield}")

        reactants = self.get_reactant_amount_pairs()
        print(f"reactants={reactants}")

        print(f"price_estimate={self.price_estimate}")

        if draw:
            self.draw()

    def draw(self) -> None:
        """Draw the molecules involved in this reaction"""

        from molparse.rdkit import draw_grid

        reactants = self.reactants

        product = self.product

        mols = [r.mol for r in reactants]
        mols.append(product.mol)

        labels = [f"+ {r}" if i > 0 else f"{r}" for i, r in enumerate(reactants)]
        labels.append(f"-> {product}")

        drawing = draw_grid(mols, labels=labels, highlightAtomLists=None)
        display(drawing)

    def check_chemistry(
        self,
        debug: bool = False,
    ) -> bool:
        """Sanity check the chemistry of this reaction

        :param debug: increase verbosity (Default value = False)
        """
        from .chem import check_chemistry

        return check_chemistry(self.type, self.reactants, self.product, debug=debug)

    def check_reactant_availability(
        self,
        supplier: None | str = None,
        debug: bool = False,
    ) -> bool:
        """Check the availability of reactant compounds

        :param supplier: Limit to quotes from this supplier (Default value = None)
        :param debug: increase verbosity (Default value = False)

        """

        if debug:
            mrich.var("reaction", self.id)
            mrich.var("reactants", self.reactant_ids)
            mrich.var("supplier", supplier)

        if supplier is None:

            triples = self.db.execute(
                f"""
                SELECT reactant_compound, SUM(quote_id), SUM(reaction_id) FROM reactant 
                LEFT JOIN quote ON quote_compound = reactant_compound
                LEFT JOIN reaction ON reaction_product = reactant_compound 
                WHERE reactant_reaction = {self.id}
                GROUP BY reactant_compound
            """
            ).fetchall()

        else:

            triples = self.db.execute(
                f"""
                WITH filtered_quotes AS
                (
                    SELECT * FROM quote
                    WHERE quote_supplier = "{supplier}"
                )
                SELECT reactant_compound, SUM(quote_id), SUM(reaction_id) FROM reactant 
                LEFT JOIN filtered_quotes ON quote_compound = reactant_compound
                LEFT JOIN reaction ON reaction_product = reactant_compound 
                WHERE reactant_reaction = {self.id}
                GROUP BY reactant_compound
            """
            ).fetchall()

        for reactant_compound, has_quote, has_reaction in triples:

            if debug:
                mrich.debug(
                    f"{reactant_compound=}, {bool(has_quote)=}, {bool(has_reaction)=}"
                )

            if has_quote:
                if debug:
                    mrich.debug(f"reactant={reactant_compound} has quote")
                continue

            if has_reaction:
                if debug:
                    mrich.debug(f"reactant={reactant_compound} has reaction")
                continue

            if debug:
                mrich.warning(f"No quote or reaction for reactant={reactant_compound}")

            return False

        return True

    def get_dict(
        self,
        smiles: bool = True,
        mols: bool = True,
    ) -> dict[str]:
        """Returns a dictionary representing this :class:`.Reaction`

        :param smiles: include smiles string (Default value = True)
        :param mols: include ``rdkit.Chem.Mol`` (Default value = True)

        """

        serialisable_fields = ["id", "type", "product_id", "reactant_ids"]

        data = {}
        for key in serialisable_fields:
            data[key] = getattr(self, key)

        if smiles:
            data["product_smiles"] = self.product_smiles
            data["reactant_smiles"] = self.reactant_smiles

        if mols:
            data["product_mol"] = self.product_mol
            data["reactant_mols"] = self.reactant_mols

        return data

    def _delete(self) -> None:
        """Delete this reaction and any related reactants, routes, and components"""

        self = animal.R53848

        route_ids = self.db.select_where(
            query="component_route",
            table="component",
            key=f"component_ref = {self.id} AND component_type = 1",
            multiple=True,
        )

        route_ids = [r for r, in route_ids]
        route_str_ids = str(tuple(route_ids)).replace(",)", ")")

        self.db.delete_where(
            table="component", key=f"component_route IN {route_str_ids}"
        )

        self.db.delete_where(table="route", key=f"route_id IN {route_str_ids}")

        self.db.delete_where(table="reactant", key="reaction", value=self.id)

        self.db.delete_where(table="reaction", key="id", value=self.id)

    ### DUNDERS

    def __str__(self) -> str:
        """Unformatted string representation"""
        return f"R{self.id}"

    def __repr__(self) -> str:
        """ANSI Formatted string representation"""
        return f"{mcol.bold}{mcol.underline}{self.plain_repr}{mcol.unbold}{mcol.ununderline}"

    def __rich__(self) -> str:
        """Rich Formatted string representation"""
        return f"[bold underline]{self.plain_repr}"

    def __eq__(
        self,
        other: "int | Reaction",
    ) -> bool:
        """compare this reaction to a :class:`.Reaction` object or ID"""

        match other:
            case int():
                return self.id == other

            case Reaction():

                if self.type != other.type:
                    return False

                if self.product != other.product:
                    return False

                if self.reactant_ids != other.reactant_ids:
                    return False

                return True

            case _:
                raise NotImplementedError

    def __hash__(self) -> int:
        """Integer hash from ID"""
        return self.id
