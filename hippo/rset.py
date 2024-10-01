from .db import Database

import mcol

import os
from numpy import int64

import logging

logger = logging.getLogger("HIPPO")

from .reaction import Reaction


class ReactionTable:
    """Class representing all :class:`.Reaction` objects in the 'reaction' table of the :class:`.Database`.

    .. attention::

            :class:`.ReactionTable` objects should not be created directly. Instead use the :meth:`.HIPPO.reactions` property. See :doc:`getting_started`.

    Use as an iterable
    ==================

    Iterate through :class:`.Reaction` objects in the table:

    ::

            for reaction in animal.reactions:
                ...


    Selecting reactions in the table
    ================================

    The :class:`.ReactionTable` can be indexed with :class:`.Reaction` ID, or list/sets/tuples/slices thereof:

    ::

            rtable = animal.reactions

            # indexing individual compounds
            reaction = rtable[13]                            # using the ID

            # getting a subset of compounds
            rset = rtable[13,15,18]      # using IDs (tuple)
            rset = rtable[[13,15,18]]    # using IDs (list)
            rset = rtable[set(13,15,18)] # using IDs (set)
            rset = rtable[13:18]         # using a slice

    """

    _name = "all reactions"

    def __init__(
        self,
        db: Database,
        table: str = "reaction",
    ) -> None:

        self._db = db
        self._table = table

    ### PROPERTIES

    @property
    def db(self) -> "Database":
        """Returns the associated :class:`.Database`"""
        return self._db

    @property
    def table(self) -> str:
        """Returns the name of the :class:`.Database` table"""
        return self._table

    @property
    def name(self) -> str | None:
        """Returns the name of set"""
        return self._name

    @property
    def types(self) -> list[str]:
        """Returns a list of the unique reaction types present in the table"""
        result = self.db.select(
            table=self.table, query="DISTINCT reaction_type", multiple=True
        )
        return [q for q, in result]

    @property
    def ids(self) -> list[int]:
        """Returns the IDs of child reactions"""
        result = self.db.select(table=self.table, query="reaction_id", multiple=True)
        return [q for q, in result]

    ### METHODS

    def interactive(self) -> None:
        """Interactive widget to navigate reactions in the table

        .. attention::

                This method instantiates a :class:`.ReactionSet` containing all poses, it is recommended to instead select a subset for display. This method is only intended for use within a Jupyter Notebook.

        """
        return self[self.ids].interactive()

    def get_by_type(self, reaction_type: str) -> "ReactionSet":
        """Get all child reactions of the given type

        :param reaction_type: reaction type to filter by

        """
        result = self.db.select_where(
            table=self.table,
            query="reaction_id",
            key="type",
            value=reaction_type,
            multiple=True,
        )
        rset = self[[q for q, in result]]
        rset._name = f"all {reaction_type} reactions"
        return rset

    def get_df(self, *, smiles: bool = True, mols: bool = True) -> "pandas.DataFrame":
        """Construct a pandas.DataFrame of all reactions in the database

        :param smiles: Include smiles column (Default value = True)
        :param mols: Include `rdkit.Chem.Mol` column (Default value = True)

        """

        from rdkit.Chem import Mol
        from pandas import DataFrame

        ### SQL QUERY

        data = {}

        if not smiles and not mols:

            sql = "SELECT reaction_id, reaction_type, reaction_product, reactant_compound FROM reaction INNER JOIN reactant ON reaction.reaction_id = reactant.reactant_reaction"

            triples = self.db.execute(sql).fetchall()

            for reaction_id, product_id, reactant_id in triples:
                if reaction_id not in data:
                    data[reaction_id] = dict(product_id=product_id, reactant_ids=[])
                else:
                    assert data[reaction_id]["product_id"] == product_id

                data[reaction_id]["reactant_ids"].append(reactant_id)

        else:

            sql = """
            SELECT {query}
            FROM reaction 

            INNER JOIN reactant 
                ON reaction.reaction_id = reactant.reactant_reaction

            INNER JOIN compound c_r
                ON c_r.compound_id = reactant.reactant_compound

            INNER JOIN compound c_p
                ON c_p.compound_id = reaction.reaction_product
            """

            if not mols:
                sql = sql.format(
                    query="reaction_id, reaction_type, reaction_product, reactant_compound, c_p.compound_smiles, c_r.compound_smiles"
                )

            else:
                sql = sql.format(
                    query="reaction_id, reaction_type, reaction_product, reactant_compound, c_p.compound_smiles, c_r.compound_smiles, mol_to_binary_mol(c_p.compound_mol), mol_to_binary_mol(c_r.compound_mol)"
                )

            results = self.db.execute(sql).fetchall()

            for result in results:

                (
                    reaction_id,
                    reaction_type,
                    product_id,
                    reactant_id,
                    product_smiles,
                    reactant_smiles,
                ) = result[:6]

                if mols:
                    product_mol, reactant_mol = result[6:]

                if reaction_id not in data:
                    data[reaction_id] = dict(
                        reaction_id=reaction_id,
                        reaction_type=reaction_type,
                        product_id=product_id,
                        reactant_ids=set(),
                        product_smiles=product_smiles,
                        reactant_smiles=set(),
                    )
                    if mols:
                        data[reaction_id]["product_mol"] = Mol(product_mol)
                        data[reaction_id]["reactant_mols"] = set()
                else:
                    assert data[reaction_id]["product_id"] == product_id

                data[reaction_id]["reactant_ids"].add(reactant_id)
                data[reaction_id]["reactant_smiles"].add(reactant_smiles)
                if mols:
                    data[reaction_id]["reactant_mols"].add(Mol(reactant_mol))

        data = data.values()
        return DataFrame(data)

    def set_product_yields(self, *, type: str, product_yield: float) -> None:
        """Set the product_yield for all member :class:`.Reaction` entries with given type

        :param type: the :class:`.Reaction` type to filter by
        :param product_yield: the :class:`.Reaction` product_yield to assign

        """

        sql = f"""
        UPDATE reaction
        SET reaction_product_yield = :reaction_product_yield
        WHERE reaction_type = :reaction_type;
        """

        self.db.execute(
            sql,
            dict(
                reaction_product_yield=product_yield,
                reaction_type=type,
            ),
        )

    ### DUNDERS

    def __getitem__(self, key) -> "Reaction | ReactionSet | None":
        """Get a member :class:`.Reaction` object or subset :class:`.ReactionSet` thereof.

        :param key: Can be an integer ID, negative integer index, list/set/tuple of IDs, or slice of IDs

        """

        match key:

            case int():

                if key == 0:
                    return self.__getitem__(key=1)

                if key < 0:
                    key = len(self) + 1 + key
                    return self.__getitem__(key=key)

                else:
                    return self.db.get_reaction(id=key)

            case key if isinstance(key, list) or isinstance(key, tuple) or isinstance(
                key, set
            ):
                return ReactionSet(self.db, key)

            case slice():
                ids = self.db.slice_ids(
                    table=self.table, start=key.start, stop=key.stop, step=key.step
                )
                return self[ids]

            case _:
                logger.error(
                    f"Unsupported type for ReactionSet.__getitem__(): {key=} {type(key)}"
                )

        return None

    def __repr__(self) -> str:
        """Formatted string representation"""

        s = f"{mcol.bold}{mcol.underline}"

        if self.name:
            s += f"{self.name}: "

        s += "{" f"R x {len(self)}" "}"

        s += f"{mcol.unbold}{mcol.ununderline}"

        return s

    def __len__(self) -> int:
        """Number of reactions in this set"""
        return self.db.count(self.table)

    def __iter__(self):
        """Iterate through poses in this set"""
        return iter(self[i + 1] for i in range(len(self)))


class ReactionSet:
    """Object representing a subset of the 'reaction' table in the :class:`.Database`.

    .. attention::

            :class:`.ReactionSet` objects should not be created directly. Instead use the :meth:`.HIPPO.reactions` property. See :doc:`getting_started` and :doc:`insert_elaborations`.

    Use as an iterable
    ==================

    Iterate through :class:`.Reaction` objects in the set:

    ::

            rset = animal.reactions[:100]

            for reaction in rset:
                    ...

    Check membership
    ================

    To determine if a :class:`.Reaction` is present in the set:

    ::

            is_member = reaction in cset

    Selecting compounds in the set
    ==============================

    The :class:`.ReactionSet` can be indexed like standard Python lists by their indices

    ::

            rset = animal.reactions[1:100]

            # indexing individual compounds
            reaction = rset[0]  # get the first reaction
            reaction = rset[1]  # get the second reaction
            reaction = rset[-1] # get the last reaction

            # getting a subset of compounds using a slice
            rset2 = rset[13:18] # using a slice

    """

    _table = "reaction"

    def __init__(
        self,
        db: Database,
        indices: list = None,
        *,
        sort: bool = True,
        name: str | None = None,
    ) -> None:

        self._db = db
        indices = indices or []

        if not isinstance(indices, list):
            indices = list(indices)

        assert all(isinstance(i, int) or isinstance(i, int64) for i in indices)

        if sort:
            self._indices = sorted(list(set(indices)))
        else:
            self._indices = list(set(indices))

        self._name = name

    ### PROPERTIES

    @property
    def db(self) -> Database:
        """Returns the associated :class:`.Database`"""
        return self._db

    @property
    def table(self) -> str:
        """Returns the name of the :class:`.Database` table"""
        return self._table

    @property
    def name(self) -> str | None:
        """Returns the name of set"""
        return self._name

    @property
    def indices(self) -> list[int]:
        """Returns the ids of poses in this set"""
        return self._indices

    @property
    def ids(self) -> list[int]:
        """Returns the ids of poses in this set"""
        return self._indices

    @property
    def types(self) -> list[str]:
        """Returns the types of poses in this set"""
        pairs = self.db.select_where(
            table="reaction",
            key=f"reaction_id IN {self.str_ids}",
            query="reaction_id, reaction_type",
            multiple=True,
        )
        lookup = {i: t for i, t in pairs}
        return [lookup[i] for i in self.ids]

    @property
    def str_ids(self) -> str:
        """Return an SQL formatted tuple string of the :class:`.Compound` IDs"""
        return str(tuple(self.ids)).replace(",)", ")")

    @property
    def products(self) -> "CompoundSet":
        """Get all product compounds that can be synthesised with these reactions (no intermediates)"""
        from .cset import CompoundSet

        intermediates = self.intermediates
        product_ids = self.db.execute(
            f"""
            SELECT compound_id FROM compound
            INNER JOIN reaction ON compound_id = reaction_product
            WHERE reaction_id IN {self.str_ids}
            AND compound_id NOT IN {intermediates.str_ids}
        """
        ).fetchall()
        cset = CompoundSet(self.db, [i for i, in product_ids])
        if self.name:
            cset._name = f"products of {self}"
        return cset

    @property
    def intermediates(self) -> "CompoundSet":
        """Get all intermediate compounds that can be synthesised with these reactions"""
        from .cset import CompoundSet

        sql = f"""
            SELECT DISTINCT compound_id FROM compound
            INNER JOIN reaction ON compound_id = reaction_product
            INNER JOIN reactant ON compound_id = reactant_compound
            WHERE reactant_reaction IN {self.str_ids}
        """
        # print(sql)
        intermediate_ids = self.db.execute(sql).fetchall()
        cset = CompoundSet(self.db, [i for i, in intermediate_ids])
        if self.name:
            cset._name = f"products of {self}"
        return cset

    ### METHODS

    def add(self, r: Reaction) -> None:
        """Add a :class:`.Reaction` to this set

        :param r: :class:`.Reaction` to be added

        """
        assert isinstance(r, Reaction)
        if (id := r.id) not in self._indices:
            self._indices.append(id)

    def interactive(self):
        """Creates a ipywidget to interactively navigate this PoseSet."""

        from ipywidgets import (
            interactive,
            BoundedIntText,
            Checkbox,
            interactive_output,
            HBox,
            GridBox,
            Layout,
            VBox,
        )
        from IPython.display import display
        from pprint import pprint

        a = BoundedIntText(
            value=0,
            min=0,
            max=len(self) - 1,
            step=1,
            description=f"Rs (/{len(self)}):",
            disabled=False,
        )

        b = Checkbox(description="Name", value=True)
        c = Checkbox(description="Summary", value=False)
        d = Checkbox(description="Draw", value=True)
        e = Checkbox(description="Check chemistry", value=False)
        f = Checkbox(description="Reactant Quotes", value=False)

        ui1 = GridBox(
            [b, c, d], layout=Layout(grid_template_columns="repeat(5, 100px)")
        )
        ui2 = GridBox([e, f], layout=Layout(grid_template_columns="repeat(2, 150px)"))
        ui = VBox([a, ui1, ui2])

        def widget(
            i, name=True, summary=True, draw=True, check_chemistry=True, reactants=False
        ):
            """

            :param i:
            :param name:  (Default value = True)
            :param summary:  (Default value = True)
            :param draw:  (Default value = True)
            :param check_chemistry:  (Default value = True)
            :param reactants:  (Default value = False)

            """
            reaction = self[i]
            if name:
                print(repr(reaction))
            if summary:
                reaction.summary(draw=False)
            if draw:
                reaction.draw()
            if check_chemistry:
                reaction.check_chemistry(debug=True)
            if reactants:
                for comp in reaction.reactants:
                    # if summary:
                    # comp.summary(draw=False)
                    # elif name:
                    print(repr(comp))

                    quotes = comp.get_quotes(df=True)
                    display(quotes)

                    # break

                    # if draw:
                    #   comp.draw()

        out = interactive_output(
            widget,
            {
                "i": a,
                "name": b,
                "summary": c,
                "draw": d,
                "check_chemistry": e,
                "reactants": f,
            },
        )

        display(ui, out)

    def get_df(self, smiles=True, mols=True, **kwargs) -> "pandas.DataFrame":
        """Construct a pandas.DataFrame of this ReactionSet

        :param smiles: Include smiles column (Default value = True)
        :param mols: Include `rdkit.Chem.Mol` column (Default value = True)
        :param kwargs: keyword arguments are passed on to :meth:`.Reaction.get_dict:

        """

        from pandas import DataFrame
        from tqdm import tqdm
        from rdkit.Chem import Mol

        logger.warning("Using slower Reaction.dict rather than direct SQL query")

        data = []
        for r in tqdm(self):
            data.append(r.get_dict(smiles=smiles, mols=mols, **kwargs))

        return DataFrame(data)

    def copy(self) -> "ReactionSet":
        """Return a copy of this set"""
        return ReactionSet(self.db, self.ids, sort=False, name=self.name)

    def get_recipes(
        self, amounts: float | list[float] = 1.0, **kwargs
    ) -> "Recipe | list[Recipe]":
        """Get the :class:`.Recipe` object(s) from this set of recipes

        :param amounts: float or list/generator of product amounts in mg, (Default value = 1.0)
        :param kwargs: keyword arguments are passed on to :meth:`.Recipe.from_reactions:

        """
        from .recipe import Recipe

        return Recipe.from_reactions(db=self.db, reactions=self, amounts=1, **kwargs)

    def reverse(self) -> None:
        """In-place reversal of indices"""
        self._indices = list(reversed(self._indices))

    def get_dict(self) -> dict[str]:
        """Serializable dictionary"""
        return dict(db=str(self.db), indices=self.indices)

    def summary(self) -> None:
        """Print a summary of the Reactions"""

        logger.header(self)
        for reaction in self:
            print(repr(reaction))

    ### DUNDERS

    def __repr__(self) -> str:
        """Formatted string representation"""

        s = f"{mcol.bold}{mcol.underline}"

        if self.name:
            s += f"{self.name}: "

        s += "{" f"R x {len(self)}" "}"

        s += f"{mcol.unbold}{mcol.ununderline}"

        return s

    def __len__(self) -> int:
        return len(self.indices)

    def __iter__(self):
        return iter(self.db.get_reaction(id=i) for i in self.indices)

    def __getitem__(self, key) -> Reaction:
        try:
            index = self.indices[key]
        except IndexError:
            logger.exception(f"list index out of range: {key=} for {self}")
            raise
        return self.db.get_reaction(id=index)

    def __add__(self, other) -> "ReactionSet":
        if other:
            for reaction in other:
                self.add(reaction)
            self._name = None
        return self
