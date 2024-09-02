import mcol

from .compound import Compound
from .recipe import Recipe

from mlog import setup_logger

logger = setup_logger("HIPPO")


class Reaction:
    """
    A :class:`.Reaction` is a simplified representation of a synthetic pathway to create a product :class:`.Compound`. Reactants (also :class:`.Compound` objects) as well as a reaction type are required.

    .. attention::

            :class:`.Reaction` objects should not be created directly. Instead use :meth:`.HIPPO.register_reaction` or :meth:`.HIPPO.reactions`

    """

    _table = "reaction"

    def __init__(
        self,
        db,
        id: int,
        type: str,
        product: int,
        product_yield: float,
    ):

        self._db = db
        self._id = id
        self._type = type
        self._product = product
        self._product_yield = product_yield

    ### FACTORIES

    ### PROPERTIES

    @property
    def id(self) -> int:
        """Returns the reaction ID"""
        return self._id

    @property
    def type(self) -> str:
        """Returns the reaction tyoe"""
        return self._type

    @property
    def product(self):
        """Returns the reaction's product :class:`.Compound`"""
        if isinstance(self._product, int):
            self._product = self.db.get_compound(id=self._product)
        return self._product

    @property
    def product_yield(self) -> float:
        """Returns the reaction's product yield (fraction)"""
        return self._product_yield

    @property
    def db(self):
        """Returns a pointer to the parent database"""
        return self._db

    @property
    def reactants(self):
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
    def reactant_ids(self) -> set:
        """Returns a set of reactant ID's"""
        return set(v for v in self.get_reactant_ids())

    @property
    def reactant_str_ids(self):
        """ """
        return str(tuple(self.reactant_ids)).replace(",)", ")")

    @property
    def product_id(self) -> int:
        """Returns the product ID"""
        return self.product.id

    @property
    def product_smiles(self):
        """ """
        return self.product.smiles

    @property
    def reactant_smiles(self):
        """ """
        return [r.smiles for r in self.reactants]

    @property
    def product_mol(self):
        """ """
        return self.product.mol

    @property
    def reactant_mols(self):
        """ """
        return [r.mol for r in self.reactants]

    @property
    def price_estimate(self):
        """ """
        return self.db.get_reaction_price_estimate(reaction=self)

    ### METHODS

    def get_reactant_amount_pairs(self, compound_object=True) -> list[Compound]:
        """Returns pairs of reactants and their amounts

        :param compound_object:  (Default value = True)

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
                    (self.db.get_compound(id=id), amount) for id, amount in compound_ids
                ]
            else:
                return compound_ids
        else:
            return []

    def get_reactant_ids(self) -> list[Compound]:
        """ """

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
        permitted_reactions=None,
        supplier: str | None = None,
    ):
        """Get a :class:`.Recipe` describing how to make the product

        :param amount: float:  (Default value = 1)
        :param debug: bool:  (Default value = False)
        :param pick_cheapest: bool:  (Default value = False)
        :param permitted_reactions:  (Default value = None)
        :param supplier: str | None:  (Default value = None)

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

    def summary(self, amount=1, draw=True):
        """Print a summary of this reaction's information

        :param amount:  (Default value = 1)
        :param draw:  (Default value = True)

        """

        print(f"id={self.id}")
        print(f"type={self.type}")
        print(f"product={self.product}")
        print(f"product_yield={self.product_yield}")

        reactants = self.get_reactant_amount_pairs()
        print(f"reactants={reactants}")

        print(f"price_estimate={self.price_estimate}")

        # print(f'Ingredients for {amount} mg of product:')
        # ingredients = self.get_ingredients(amount=amount)
        # print(ingredients)

        if draw:
            self.draw()

        # return self.get_recipe(amount)

    def draw(self):
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

    def check_chemistry(self, debug=False):
        """

        :param debug:  (Default value = False)

        """
        from .chem import check_chemistry

        return check_chemistry(self.type, self.reactants, self.product, debug=debug)

    def check_reactant_availability(
        self, supplier: None | str = None, debug: bool = False
    ):
        """

        :param supplier: None | str:  (Default value = None)
        :param debug: bool:  (Default value = False)

        """

        if debug:
            logger.var("reaction", self.id)
            logger.var("reactants", self.reactant_ids)
            logger.var("supplier", supplier)

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
                logger.debug(
                    f"{reactant_compound=}, {bool(has_quote)=}, {bool(has_reaction)=}"
                )

            if has_quote:
                if debug:
                    logger.debug(f"reactant={reactant_compound} has quote")
                continue

            if has_reaction:
                if debug:
                    logger.debug(f"reactant={reactant_compound} has reaction")
                continue

            if debug:
                logger.warning(f"No quote or reaction for reactant={reactant_compound}")

            return False

        return True

    def get_dict(self, smiles=True, mols=True):
        """Returns a dictionary representing this Reaction

        :param smiles:  (Default value = True)
        :param mols:  (Default value = True)

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

    ### DUNDERS

    def __str__(self):
        return f"R{self.id}"

    def __repr__(self):
        return f"{mcol.bold}{mcol.underline}{self}: {self.reaction_str} via {self.type}{mcol.unbold}{mcol.ununderline}"

    def __eq__(self, other):

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

    def __hash__(self):
        return self.id
