from .db import Database

import mcol

import os
from numpy import int64

import logging

logger = logging.getLogger("HIPPO")

from .reaction import Reaction


class ReactionTable:
    """Object representing the 'reaction' table in the :class:`.Database`."""

    def __init__(
        self,
        db: Database,
        table: str = "reaction",
    ) -> None:

        self._db = db
        self._table = table

    ### PROPERTIES

    @property
    def db(self):
        """Returns the associated :class:`.Database`"""
        return self._db

    @property
    def table(self):
        """ """
        return self._table

    @property
    def types(self):
        """ """
        result = self.db.select(
            table=self.table, query="DISTINCT reaction_type", multiple=True
        )
        return [q for q, in result]

    @property
    def ids(self):
        """Returns the IDs of child reactions"""
        result = self.db.select(table=self.table, query="reaction_id", multiple=True)
        return [q for q, in result]

    ### METHODS

    def interactive(self):
        """ """
        return self[self.ids].interactive()

    def get_by_type(self, reaction_type: str):
        """

        :param reaction_type: str:

        """
        result = self.db.select_where(
            table=self.table,
            query="reaction_id",
            key="type",
            value=reaction_type,
            multiple=True,
        )
        return self[[q for q, in result]]

    def get_df(self, smiles=True, mols=True, **kwargs):
        """Construct a pandas.DataFrame of all reactions in the database

        :param smiles:  (Default value = True)
        :param mols:  (Default value = True)

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

    ### DUNDERS

    def __getitem__(self, key) -> Reaction:

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
        return f"{mcol.bold}{mcol.underline}set(R x {len(self)}){mcol.unbold}{mcol.ununderline}"

    def __len__(self) -> int:
        return self.db.count(self.table)

    def __iter__(self):
        return iter(self[i + 1] for i in range(len(self)))


class ReactionSet:
    """Object representing a subset of the 'reaction' table in the :class:`.Database`."""

    _table = "reaction"

    def __init__(
        self,
        db: Database,
        indices: list = None,
        *,
        sort: bool = True,
    ):

        self._db = db
        indices = indices or []

        if not isinstance(indices, list):
            indices = list(indices)

        assert all(isinstance(i, int) or isinstance(i, int64) for i in indices)

        if sort:
            self._indices = sorted(list(set(indices)))
        else:
            self._indices = list(set(indices))

    ### PROPERTIES

    @property
    def db(self):
        """Returns the associated :class:`.Database`"""
        return self._db

    @property
    def table(self):
        """ """
        return self._table

    @property
    def indices(self):
        """ """
        return self._indices

    @property
    def ids(self):
        """ """
        return self._indices

    @property
    def str_ids(self):
        """ """
        return str(tuple(self.ids)).replace(",)", ")")

    @property
    def products(self):
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
        return CompoundSet(self.db, [i for i, in product_ids])

    @property
    def intermediates(self):
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
        return CompoundSet(self.db, [i for i, in intermediate_ids])

    ### METHODS

    def add(self, r):
        """

        :param r:

        """
        assert r._table == "reaction"
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

    def get_df(self, smiles=True, mols=True, **kwargs):
        """Construct a pandas.DataFrame of this ReactionSet

        :param smiles:  (Default value = True)
        :param mols:  (Default value = True)

        """

        from pandas import DataFrame
        from tqdm import tqdm
        from rdkit.Chem import Mol

        logger.warning("Using slower Reaction.dict rather than direct SQL query")

        data = []
        for r in tqdm(self):
            data.append(r.get_dict(smiles=smiles, mols=mols, **kwargs))

        return DataFrame(data)

    def copy(self):
        """ """
        return ReactionSet(self.db, self.ids, sort=False)

    def get_recipes(self, amounts=1):
        """

        :param amounts:  (Default value = 1)

        """
        from .recipe import Recipe

        return Recipe.from_reactions(db=self.db, reactions=self, amounts=1)

    def reverse(self):
        """ """
        self._indices = list(reversed(self._indices))

    def get_dict(self):
        """Serializable dictionary"""
        return dict(db=str(self.db), indices=self.indices)

    ### DUNDERS

    def __repr__(self) -> str:
        return (
            f"{mcol.bold}{mcol.underline}"
            "{"
            f"R x {len(self)}"
            "}"
            f"{mcol.unbold}{mcol.ununderline}"
        )

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

    def __add__(self, other):
        for reaction in other:
            self.add(reaction)
        return self
