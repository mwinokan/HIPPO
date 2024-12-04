# from .tools import df_row_to_dict

from .compound import Compound, Ingredient
from .db import Database

from .recipe import Recipe

from numpy import int64, nan, isnan, mean, std

import mcol

import os

import mrich


class CompoundTable:
    """Class representing all :class:`.Compound` objects in the 'compound' table of the :class:`.Database`.

    .. attention::

            :class:`.CompoundTable` objects should not be created directly. Instead use the :meth:`.HIPPO.compounds` property. See :doc:`getting_started` and :doc:`insert_elaborations`.

    Use as an iterable
    ==================

    Iterate through :class:`.Compound` objects in the table:

    ::

            for compound in animal.compounds:
                    ...


    Selecting compounds in the table
    ================================

    The :class:`.CompoundTable` can be indexed with :class:`.Compound` IDs, names, aliases, or list/sets/tuples/slices thereof:

    ::

            ctable = animal.compounds

            # indexing individual compounds
            comp = ctable[13]                            # using the ID
            comp = ctable["BSYNRYMUTXBXSQ-UHFFFAOYSA-N"] # using the InChIKey
            comp = ctable["aspirin"]                     # using the alias

            # getting a subset of compounds
            cset = ctable[13,15,18]      # using IDs (tuple)
            cset = ctable[[13,15,18]]    # using IDs (list)
            cset = ctable[set(13,15,18)] # using IDs (set)
            cset = ctable[13:18]         # using a slice

    Tags and base compounds can also be used to filter:

    ::

            cset = animal.compounds(tag='hits') # select compounds tagged with 'hits'
            cset = animal.compounds(base=comp)  # select elaborations of comp

    """

    _table = "compound"
    _name = "all compounds"

    def __init__(
        self,
        db: Database,
    ) -> None:

        self._db = db

    ### FACTORIES

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
    def names(self) -> list[str]:
        """Returns the names of child compounds"""
        result = self.db.select(table=self.table, query="compound_name", multiple=True)
        return [q for q, in result]

    @property
    def name(self):
        return self._name

    @property
    def ids(self) -> list[int]:
        """Returns the IDs of child compounds"""
        result = self.db.select(table=self.table, query="compound_id", multiple=True)
        return [q for q, in result]

    @property
    def tags(self) -> set[str]:
        """Returns the set of unique tags present in this compound set"""
        values = self.db.select_where(
            table="tag",
            query="DISTINCT tag_name",
            key="tag_compound IS NOT NULL",
            multiple=True,
        )
        return set(v for v, in values)

    @property
    def reactants(self) -> "CompoundSet":
        """Returns a :class:`.CompoundSet` of all compounds that are used as a reactants"""
        # ids = self.db.select(table='reactant', query='DISTINCT reactant_compound', multiple=True)
        ids = self.db.execute(
            "SELECT reactant_compound FROM reactant LEFT JOIN reaction ON reactant.reactant_compound = reaction.reaction_product WHERE reaction.reaction_product IS NULL"
        ).fetchall()
        ids = [q for q, in ids]
        from .cset import CompoundSet

        cset = CompoundSet(self.db, ids)
        cset._name = "all reactants"
        return cset

    @property
    def products(self) -> "CompoundSet":
        """Returns a :class:`.CompoundSet` of all compounds that are a product of a reaction but not a reactant"""
        ids = self.db.execute(
            "SELECT reaction_product FROM reaction LEFT JOIN reactant ON reaction.reaction_product = reactant.reactant_compound WHERE reactant.reactant_compound IS NULL"
        ).fetchall()
        ids = [q for q, in ids]
        from .cset import CompoundSet

        cset = CompoundSet(self.db, ids)
        cset._name = "all products"
        return cset

    @property
    def intermediates(self) -> "CompoundSet":
        """Returns a :class:`.CompoundSet` of all compounds that are products and reactants"""
        ids = self.db.execute(
            "SELECT DISTINCT reaction_product FROM reaction INNER JOIN reactant ON reaction.reaction_product = reactant.reactant_compound"
        ).fetchall()
        ids = [q for q, in ids]
        from .cset import CompoundSet

        cset = CompoundSet(self.db, ids)
        cset._name = "all intermediates"
        return cset

    @property
    def num_reactants(self) -> int:
        """Returns the number of reactants (see :meth:`CompoundTable.reactants`)"""
        return len(self.reactants)

    @property
    def num_intermediates(self) -> int:
        """Returns the number of intermediates (see :meth:`CompoundTable.intermediates`)"""
        return len(self.intermediates)

    @property
    def num_products(self) -> int:
        """Returns the number of products (see :meth:`CompoundTable.products`)"""
        return len(self.products)

    @property
    def elabs(self) -> "CompoundSet":
        """Returns a :class:`.CompoundSet` of all compounds that are a an elaboration of an existing base"""
        ids = self.db.select_where(
            query="scaffold_superstructure",
            table="scaffold",
            key="scaffold_superstructure IS NOT NULL",
            multiple=True,
            none="quiet",
        )

        if not ids:
            return None

        ids = [q for q, in ids]
        from .cset import CompoundSet

        cset = CompoundSet(self.db, ids)
        cset._name = "all elaborations"
        return cset

    @property
    def bases(self) -> "CompoundSet":
        """Returns a :class:`.CompoundSet` of all compounds that are the basis for a set of elaborations"""
        ids = self.db.select_where(
            query="DISTINCT scaffold_base",
            table="scaffold",
            key="scaffold_base IS NOT NULL",
            multiple=True,
            none="quiet",
        )
        ids = [q for q, in ids]
        from .cset import CompoundSet

        cset = CompoundSet(self.db, ids)
        cset._name = "all bases"
        return cset

    @property
    def num_elabs(self) -> int:
        """Returns the number of compounds that are a an elaboration of an existing base"""
        return len(self.elabs)

    @property
    def num_bases(self) -> int:
        """Returns the number of compounds that are the basis for a set of elaborations"""
        return len(self.bases)

    ### METHODS

    def get_by_tag(
        self,
        tag: str,
        **kwargs,
    ) -> "CompoundSet":
        """Get all child compounds with a certain tag

        :param tag: tag to filter by

        """
        values = self.db.select_where(
            query="tag_compound",
            table="tag",
            key="name",
            value=tag,
            multiple=True,
            **kwargs,
        )

        if not values:
            return None

        ids = [v for v, in values if v]
        cset = self[ids]
        cset._name = f"compounds tagged {tag}"
        return cset

    def get_by_metadata(
        self,
        key: str,
        value: str | None = None,
    ):
        """Get all child compounds by their metadata. If no value is passed, then simply containing the key in the metadata dictionary is sufficient

        :param key: metadata key
        :param value: metadata value (Default value = None)

        """
        results = self.db.select(
            query="compound_id, compound_metadata", table="compound", multiple=True
        )
        if value is None:
            ids = [i for i, d in results if d and f'"{key}":' in d]
            name = f"compounds with {key} in metadata"
        else:
            if isinstance(value, str):
                value = f'"{value}"'
            ids = [i for i, d in results if d and f'"{key}": {value}' in d]
            name = f"compounds with metadata[{key}] == {value}"

        cset = self[ids]
        cset._name = name
        return cset

    def get_by_base(
        self,
        base: Compound | int,
    ) -> "CompoundSet":
        """Get all compounds that elaborate the given base compound

        :param base: :class:`.Compound` object or ID to search by

        """

        if not isinstance(base, int):
            assert base._table == "compound"
            base = base.id

        values = self.db.select_where(
            query="scaffold_superstructure",
            table="scaffold",
            key="base",
            value=base,
            multiple=True,
        )
        ids = [v for v, in values if v]
        cset = self[ids]
        cset._name = f"elaborations of C{base}"
        return cset

    def get_by_smiles(self, smiles: str, **kwargs) -> "Compound | None":
        """Get a member compound by its smiles"""

        from .tools import inchikey_from_smiles, sanitise_smiles, SanitisationError

        assert isinstance(smiles, str), f"Non-string {smiles=}"
        try:
            smiles = sanitise_smiles(smiles, sanitisation_failed="error")
        except SanitisationError as e:
            mrich.error(f"Could not sanitise {smiles=}")
            mrich.error(str(e))
            return None
        except AssertionError:
            mrich.error(f"Could not sanitise {smiles=}")
            return None
            return c
        inchikey = inchikey_from_smiles(smiles)
        return self.db.get_compound(inchikey=inchikey, **kwargs)

    def summary(self) -> None:
        """Print a summary of this compound set"""
        mrich.header("CompoundTable()")
        mrich.var("#compounds", len(self))
        # mrich.var('#poses', self.num_poses)
        mrich.var("tags", self.tags)
        mrich.var("#bases", self.num_bases)
        mrich.var("#elabs", self.num_elabs)
        mrich.var("#reactants", self.num_reactants)
        mrich.var("#intermediates", self.num_intermediates)
        mrich.var("#products", self.num_products)

    def draw(self) -> None:
        """2D grid of drawings of molecules in this set

        .. attention::

                This method instantiates a :class:`.CompoundSet` containing all compounds, it is recommended to instead select a subset for display. This method is only intended for use within a Jupyter Notebook.

        """
        return self[self.ids].draw()

    def interactive(self) -> None:
        """Interactive widget to navigate compounds in the table

        .. attention::

                This method instantiates a :class:`.CompoundSet` containing all compounds, it is recommended to instead select a subset for display. This method is only intended for use within a Jupyter Notebook.

        """
        self[self.ids].interactive()

    ### DUNDERS

    def __call__(
        self,
        *,
        tag: str = None,
        base: int | Compound = None,
        smiles: str | None = None,
        **kwargs,
    ) -> "CompoundSet | Compound | None":
        """Filter compounds by a given tag, base, or it's SMILES string. See :meth:`.CompoundTable.get_by_tag` and :meth:`.CompoundTable.get_by_base`

        :param tag: optional tag to filter by
        :param base: optional :class:`.Compound` ID or object to filter by
        :param base: optional SMILES string to filter by
        :returns: :class:`.CompoundSet` if searching by tag or base, else :class:`.Compound` object

        """

        if tag:
            return self.get_by_tag(tag, **kwargs)
        elif base:
            return self.get_by_base(base, **kwargs)
        elif smiles:
            return self.get_by_smiles(smiles, **kwargs)
        else:
            mrich.error("Must provide one of tag, base, or smiles arguments")
            return None

    def __getitem__(
        self,
        key: int | str | tuple | list | set | slice,
    ) -> Compound:
        """Get a member :class:`.Pose` object or subset :class:`.PoseSet` thereof.

        :param key: Can be an integer ID, negative integer index, alias or inchikey string, list/set/tuple of IDs, or slice of IDs

        """

        match key:

            # case int():
            case key if isinstance(key, int) or isinstance(key, int64):

                if key == 0:
                    return self.__getitem__(key=1)

                if key < 0:
                    key = len(self) + 1 + key
                    return self.__getitem__(key=key)

                else:
                    return self.db.get_compound(id=key)

            case str():
                comp = self.db.get_compound(inchikey=key)
                if not comp:
                    comp = self.db.get_compound(alias=key)
                return comp

            case key if isinstance(key, list) or isinstance(key, tuple) or isinstance(
                key, set
            ):

                indices = []
                for i in key:
                    if isinstance(i, int) or isinstance(i, int64):
                        index = i
                    elif isinstance(i, str):
                        index = self.db.get_compound_id(name=i)
                    else:
                        raise NotImplementedError

                    assert index
                    indices.append(index)

                return CompoundSet(self.db, indices)

            case slice():
                ids = self.db.slice_ids(
                    table=self.table, start=key.start, stop=key.stop, step=key.step
                )
                return self[ids]

            case _:
                mrich.error(
                    f"Unsupported type for CompoundTable.__getitem__(): {key=} {type(key)}"
                )

        return None

    def __str__(self) -> str:
        """Unformatted string representation"""

        if self.name:
            s = f"{self.name}: "
        else:
            s = ""

        s += "{" f"C × {len(self)}" "}"

        return s

    def __repr__(self) -> str:
        """ANSI ormatted string representation"""
        return f"{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}"

    def __rich__(self) -> str:
        """Representation for mrich"""
        return f"[bold underline]{self}"

    def __len__(self) -> int:
        """Total number of compounds"""
        return self.db.count(self.table)

    def __iter__(self):
        """Iterate through all compounds"""
        return iter(self[i + 1] for i in range(len(self)))


class CompoundSet:
    """Object representing a subset of the 'compound' table in the :class:`.Database`.

    .. attention::

            :class:`.CompoundSet` objects should not be created directly. Instead use the :meth:`.HIPPO.compounds` property. See :doc:`getting_started` and :doc:`insert_elaborations`.

    Use as an iterable
    ==================

    Iterate through :class:`.Compound` objects in the set:

    ::

            cset = animal.compounds[:100]

            for compound in cset:
                    ...

    Check membership
    ================

    To determine if a :class:`.Compound` is present in the set:

    ::

            is_member = compound in cset

    Selecting compounds in the set
    ==============================

    The :class:`.CompoundSet` can be indexed like standard Python lists by their indices

    ::

            cset = animal.compounds[1:100]

            # indexing individual compounds
            comp = cset[0]  # get the first compound
            comp = cset[1]  # get the second compound
            comp = cset[-1] # get the last compound

            # getting a subset of compounds using a slice
            cset2 = cset[13:18] # using a slice

    Tags and base compounds can also be used to filter:

    ::

            cset = animal.compounds(tag='hits') # select compounds tagged with 'hits'
            cset = animal.compounds(base=comp)  # select elaborations of comp

    """

    _table = "compound"

    def __init__(
        self,
        db: Database,
        indices: list = None,
        sort: bool = True,
        name: str | None = None,
    ):

        self._db = db

        indices = indices or []

        if not isinstance(indices, list):
            indices = list(indices)

        indices = [int(i) for i in indices]

        if sort:
            self._indices = sorted(list(set(indices)))
        else:
            self._indices = list(set(indices))

        self._name = name
        self._total_changes = db.total_changes

    ### PROPERTIES

    @property
    def db(self) -> "Database":
        """ """
        return self._db

    @property
    def table(self) -> str:
        """Get the name of the database table"""
        return self._table

    @property
    def indices(self) -> list[int]:
        """Returns the ids of compounds in this set"""
        return self._indices

    @property
    def ids(self) -> list[int]:
        """Returns the ids of compounds in this set"""
        return self.indices

    @property
    def name(self) -> str | None:
        """Returns the name of set"""
        return self._name

    @property
    def names(self) -> list[str]:
        """Returns the aliases of compounds in this set"""
        result = self.db.select_where(
            query="compound_alias",
            table="compound",
            key=f"compound_id in {self.str_ids}",
            multiple=True,
        )
        return [q for q, in result]

    @property
    def smiles(self) -> list[str]:
        """Returns the smiles of child compounds"""
        result = self.db.select_where(
            query="compound_smiles",
            table="compound",
            key=f"compound_id in {self.str_ids}",
            multiple=True,
        )
        return [q for q, in result]

    @property
    def inchikeys(self) -> list[str]:
        """Returns the inchikeys of compounds in this set"""
        result = self.db.select_where(
            query="compound_inchikey",
            table="compound",
            key=f"compound_id in {self.str_ids}",
            multiple=True,
        )
        return [q for q, in result]

    @property
    def tags(self) -> set[str]:
        """Returns the set of unique tags present in this compound set"""
        values = self.db.select_where(
            table="tag",
            query="DISTINCT tag_name",
            key=f"tag_compound in {self.str_ids}",
            multiple=True,
        )
        if not values:
            return set()
        return set(v for v, in values)

    @property
    def num_poses(self) -> int:
        """Count the poses associated to this set of compounds"""
        from .pset import PoseSet

        return self.db.count_where(table="pose", key=f"pose_compound in {self.str_ids}")

    @property
    def poses(self) -> "PoseSet":
        """Get the poses associated to this set of compounds"""
        from .pset import PoseSet

        ids = self.db.select_where(
            query="pose_id",
            table="pose",
            key=f"pose_compound in {self.str_ids}",
            multiple=True,
        )
        ids = [v for v, in ids]
        return PoseSet(self.db, ids)

    @property
    def best_placed_poses(self) -> "PoseSet":
        """Get the best placed pose for each compound in this set"""
        from .pset import PoseSet

        query = self.db.select_where(
            table="pose",
            query="pose_id, MIN(pose_distance_score)",
            key=f"pose_compound in {self.str_ids} GROUP BY pose_compound",
            multiple=True,
        )
        ids = [i for i, s in query]
        return PoseSet(self.db, ids)

    @property
    def str_ids(self) -> str:
        """Return an SQL formatted tuple string of the :class:`.Compound` IDs"""
        return str(tuple(self.ids)).replace(",)", ")")

    @property
    def num_heavy_atoms(self) -> int:
        """Get the total number of heavy atoms"""
        return sum([c.num_heavy_atoms for c in self])

    @property
    def num_rings(self):
        """Get the total number of molecular rings"""
        return sum([c.num_rings for c in self])

    @property
    def formula(self) -> str:
        """Get the combined chemical formula for all compounds"""
        from molparse.atomtypes import atomtype_dict_to_formula

        return atomtype_dict_to_formula(self.atomtype_dict)

    @property
    def atomtype_dict(self) -> dict[str, int]:
        """Get a dictionary with atomtypes as keys and corresponding quantities/counts as values"""
        from molparse.atomtypes import formula_to_atomtype_dict, combine_atomtype_dicts

        atomtype_dicts = [c.atomtype_dict for c in self]
        return combine_atomtype_dicts(atomtype_dicts)

    @property
    def num_atoms_added(self) -> list[int]:
        """Calculate the number of atoms added w.r.t the base

        :returns: list of number of atoms added values

        """

        query = self.db.execute(
            f"""
        WITH nums AS (
            SELECT A.compound_id AS comp_id, 
            mol_num_hvyatms(A.compound_mol) - mol_num_hvyatms(B.compound_mol) AS diff 
            FROM compound A, compound B
            WHERE A.compound_base = B.compound_id
            AND A.compound_id IN {self.str_ids}
        )

        SELECT compound_id, diff FROM compound
        LEFT JOIN nums
        ON comp_id = compound_id
        WHERE compound_id IN {self.str_ids}
        """
        ).fetchall()

        lookup = {k: v for k, v in query}

        return [lookup[i] for i in self.indices]

    @property
    def avg_num_atoms_added(self) -> float:
        """Calculate the average number of atoms added w.r.t the base

        :returns: average number of atoms added values for compounds which have a base

        """

        (avg,) = self.db.execute(
            f"""
        WITH nums AS (
            SELECT A.compound_id AS comp_id, 
            mol_num_hvyatms(A.compound_mol) - mol_num_hvyatms(B.compound_mol) AS diff 
            FROM compound A, compound B
            WHERE A.compound_base = B.compound_id
            AND A.compound_id IN {self.str_ids}
        )

        SELECT AVG(diff) FROM compound
        INNER JOIN nums
        ON comp_id = compound_id
        WHERE compound_id IN {self.str_ids}
        """
        ).fetchone()

        return avg

    @property
    def risk_diversity(self) -> float:
        """Calculate the average spread of risk (#atoms added) for each base in this set

        :returns: average of the standard deviations of number of atoms added for each base

        """

        return self.get_risk_diversity()

    @property
    def elaboration_balance(self) -> float:
        """Measure of how evenly elaborations are distributed across bases in this set"""

        sql = f"""
        SELECT COUNT(1) FROM scaffold
        WHERE scaffold_superstructure IN {self.str_ids}
        GROUP BY scaffold_base
        """

        counts = self.db.execute(sql).fetchall()

        counts = [c for c, in counts]  # + [0 for _ in range(len(self)-len(counts))]

        from hirsch import hirsch

        return hirsch(counts)

        # return -std(counts)

    @property
    def num_bases_elaborated(self) -> int:
        """Count the number of base compounds that have at least one elaboration in this set

        :returns: number of base compounds

        """

        (count,) = self.db.execute(
            f"""
                SELECT COUNT(DISTINCT scaffold_base) FROM scaffold
                WHERE scaffold_superstructure IN {self.str_ids}  
            """
        ).fetchone()

        return count

    @property
    def bases(self) -> "CompoundSet":
        """Get the base compounds that have at least one elaboration in this set

        :returns: :class:`.CompoundSet`

        """

        base_ids = self.db.execute(
            f"""
                SELECT DISTINCT scaffold_base FROM scaffold
                WHERE scaffold_superstructure IN {self.str_ids}  
            """
        ).fetchall()

        base_ids = [i for i, in base_ids]
        return CompoundSet(self.db, base_ids)

    @property
    def elabs(self) -> "CompoundSet":
        """Returns a :class:`.CompoundSet` of all compounds that are a an elaboration of an existing base"""

        ids = self.db.select_where(
            query="scaffold_superstructure",
            table="scaffold",
            key=f"scaffold_superstructure IS NOT NULL and scaffold_base IN {self.str_ids}",
            multiple=True,
            none="quiet",
        )

        if not ids:
            return None

        ids = [q for q, in ids]
        from .cset import CompoundSet

        return CompoundSet(self.db, ids)

    @property
    def id_num_poses_dict(self) -> dict[int, int]:
        """Get a dictionary mapping compound ids to the number of poses"""

        sql = f"""
            SELECT pose_compound, COUNT(1) FROM pose
            WHERE pose_compound IN {self.str_ids}
            GROUP BY pose_compound
        """

        records = self.db.execute(sql)

        assert records

        lookup = {k: v for k, v in records}

        for id in self.ids:
            if id not in lookup:
                lookup[id] = 0

        return lookup

    @property
    def _db_changed(self) -> bool:
        """Has the database changed?"""
        if self._total_changes != self.db.total_changes:
            self._total_changes = self.db.total_changes
            return True
        return False

    @property
    def reaction_ids(self) -> list[int]:
        """Returns a list of :class:`.Reaction` IDs that result in members of this set"""
        records = self.db.select_where(
            table="reaction",
            query="reaction_id",
            key=f"reaction_product IN {self.str_ids}",
            multiple=True,
        )
        if not records:
            return None
        return [r for r, in records]

    ### FILTERING

    def get_by_tag(
        self,
        tag: str,
    ) -> "CompoundSet":
        """Get all child compounds with a certain tag"""

        values = self.db.select_where(
            query="tag_compound",
            table="tag",
            key=f'tag_name = "{tag}" AND tag_compound IN {self.str_ids}',
            multiple=True,
        )
        ids = [v for v, in values if v and v in self.ids]
        return CompoundSet(self.db, ids)

    def get_by_metadata(self, key: str, value: str | None = None) -> "CompoundSet":
        """Get all child compounds with by their metadata. If no value is passed, then simply containing the key in the metadata dictionary is sufficient

        :param key: metadata key
        :param value: metadata value (Default value = None)
        """

        results = self.db.select(
            query="compound_id, compound_metadata", table="compound", multiple=True
        )
        if value is None:
            ids = [i for i, d in results if d and f'"{key}":' in d and i in self.ids]
        else:
            if isinstance(value, str):
                value = f'"{value}"'
            ids = [
                i
                for i, d in results
                if d and f'"{key}": {value}' in d and i in self.ids
            ]
        return CompoundSet(self.db, ids)

    def get_by_base(
        self,
        base: Compound | int,
        none: str = "error",
    ) -> "CompoundSet":
        """Get all compounds that elaborate the given base compound

        :param base: :class:`.Compound` object or ID to search by

        """

        if not isinstance(base, int):
            assert base._table == "compound"
            base = base.id

        values = self.db.select_where(
            query="scaffold_superstructure",
            table="scaffold",
            key=f"scaffold_base = {base} AND scaffold_superstructure IN {self.str_ids}",
            multiple=True,
            none=none,
        )
        ids = [v for v, in values if v]

        if not ids:
            return None
        return CompoundSet(self.db, ids)

    def get_all_possible_reactants(
        self,
        debug: bool = False,
    ) -> "CompoundSet":
        """Recursively searches for all the reactants that could possible be needed to synthesise these compounds.

        :param debug: Increased verbosity for debugging (Default value = False)

        """
        all_reactants, all_reactions = self.db.get_unsolved_reaction_tree(
            product_ids=self.ids, debug=debug
        )
        return all_reactants

    def get_all_possible_reactions(
        self,
        debug: bool = False,
    ) -> "ReactionSet":
        """Recursively searches for all the reactants that could possible be needed to synthesise these compounds.

        :param debug: Increased verbosity for debugging (Default value = False)

        """
        all_reactants, all_reactions = self.db.get_unsolved_reaction_tree(
            product_ids=self.ids, debug=debug
        )
        return all_reactions

    def get_risk_diversity(self, debug: bool = False) -> float:
        """Calculate the average spread of risk (#atoms added) for each base in this set

        :returns: average of the standard deviations of number of atoms added for each base

        """

        variances = self.db.execute(
            f"""
        WITH nums AS (
            SELECT scaffold_base AS base, scaffold_superstructure AS elab, 
            mol_num_hvyatms(c2.compound_mol) - mol_num_hvyatms(c1.compound_mol) AS diff
            FROM scaffold
            INNER JOIN compound AS c1 ON scaffold_base = c1.compound_id
            INNER JOIN compound AS c2 ON scaffold_superstructure = c2.compound_id
            WHERE scaffold_superstructure IN {self.str_ids}
        ),

        means AS (  
            SELECT base, AVG(diff) AS mean FROM nums
            GROUP BY base
        )

        SELECT AVG((nums.diff - mean)*(nums.diff - mean)) var FROM nums
        LEFT JOIN means
        ON nums.base = means.base
        GROUP BY nums.base
        """
        ).fetchall()

        if not variances:
            return None

        variances = [v for v, in variances]

        if debug:
            mrich.debug(f"{variances=}")

        return mean(variances)

    def count_by_tag(
        self,
        tag: str,
    ) -> "CompoundSet":
        """Count all child compounds with a certain tag

        :param tag: tag to filter by

        """
        (count,) = self.db.select_where(
            query="COUNT(tag_compound)",
            table="tag",
            key=f'tag_name = "{tag}" AND tag_compound IN {self.str_ids}',
            multiple=False,
        )
        return count

    ### CONSOLE / NOTEBOOK OUTPUT

    def draw(self) -> None:
        """Draw a grid of all contained molecules.

        .. attention::

                This method is only intended for use within a Jupyter Notebook.

        """

        from molparse.rdkit import draw_grid

        data = [(str(c), c.mol) for c in self]

        mols = [d[1] for d in data]
        labels = [d[0] for d in data]

        display(draw_grid(mols, labels=labels))

    def grid(self) -> None:
        """Draw a grid of all contained molecules.

        .. attention::

                This method is only intended for use within a Jupyter Notebook.

        """

        self.draw()

    def summary(self) -> None:
        """Print a summary of this compound set"""
        mrich.header("CompoundSet()")
        mrich.var("#compounds", len(self))
        mrich.var("#poses", self.num_poses)
        mrich.var("tags", self.tags)

    def interactive(self) -> None:
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
            description=f"Comp (/{len(self)}):",
            disabled=False,
        )

        b = Checkbox(description="Name", value=True)
        c = Checkbox(description="Summary", value=False)
        d = Checkbox(description="2D", value=True)
        e = Checkbox(description="poses", value=False)
        f = Checkbox(description="reactions", value=False)
        g = Checkbox(description="Metadata", value=False)

        ui = GridBox(
            [b, c, d, e, f, g], layout=Layout(grid_template_columns="repeat(6, 100px)")
        )
        ui = VBox([a, ui])

        def widget(
            i,
            name=True,
            summary=True,
            draw=True,
            poses=True,
            reactions=True,
            metadata=True,
        ):
            """

            :param i: param name:  (Default value = True)
            :param summary: Default value = True)
            :param draw: Default value = True)
            :param poses: Default value = True)
            :param reactions: Default value = True)
            :param metadata: Default value = True)
            :param name:  (Default value = True)

            """
            comp = self[i]

            if name and not summary:
                print(repr(comp))

            if summary:
                comp.summary(metadata=False, draw=False)

            if draw:
                comp.draw()

            if poses and (pset := comp.poses):
                for p in pset:
                    print(repr(p))
                pset.draw()

            if reactions and (reactions := comp.reactions):
                for r in reactions:
                    print(repr(r))
                    r.draw()

            if metadata:
                mrich.title("Metadata:")
                pprint(comp.metadata)

        out = interactive_output(
            widget,
            {
                "i": a,
                "name": b,
                "summary": c,
                "draw": d,
                "poses": e,
                "reactions": f,
                "metadata": g,
            },
        )

        display(ui, out)

    ### OTHER METHODS

    def add(self, compound: Compound | int) -> None:
        """Add a compound to this set

        :param compound: compound to be added

        """

        if isinstance(compound, Compound):
            compound = compound.id

        if compound not in self.ids:
            from bisect import insort

            insort(self.ids, compound)

    def get_recipes(
        self,
        amount: float = 1,
        debug: bool = False,
        pick_cheapest: bool = False,
        permitted_reactions: "ReactionSet | None" = None,
        quoted_only: bool = False,
        supplier: None | str = None,
        **kwargs,
    ):
        """Generate the :class:`.Recipe` to make these compounds.

        See :meth:`.Recipe.from_compounds`
        """
        from .recipe import Recipe

        return Recipe.from_compounds(
            self,
            amount=amount,
            debug=debug,
            pick_cheapest=pick_cheapest,
            permitted_reactions=permitted_reactions,
            quoted_only=quoted_only,
            supplier=supplier,
            **kwargs,
        )

    def get_routes(
        self,
        permitted_reactions: "ReactionSet",
        debug: bool = True,
    ):

        if "route" not in self.db.table_names:
            mrich.error("route table not in Database")
            raise NotImplementedError

        sql = f"""
        SELECT route_id, route_product, component_ref FROM route
        INNER JOIN component ON route_id = component_route
        WHERE route_product IN {self.str_ids}
        AND component_type = 1
        """

        permitted_reactions = set(permitted_reactions.ids)

        if debug:
            mrich.debug("Querying database for routes")
        records = self.db.execute(sql).fetchall()

        if debug:
            mrich.debug("Assembling route dictionary")

        routes = {}
        for route_id, route_product, reaction_id in records:
            if route_id not in routes:
                routes[route_id] = dict(product=route_product, reactions=set())
            assert routes[route_id]["product"] == route_product
            routes[route_id]["reactions"].add(reaction_id)

        if debug:
            mrich.debug("Checking availability")
        available_routes = set()
        for route_id, route_dict in routes.items():
            product = route_dict["product"]
            assert product in self
            reactions = route_dict["reactions"]
            if all(r in permitted_reactions for r in reactions):
                available_routes.add(route_id)

        routes = [
            self.db.get_route(id=route_id)
            for route_id in mrich.track(available_routes, prefix="Getting routes")
        ]

        from .recipe import RouteSet

        return RouteSet(self.db, routes)

    def copy(self) -> "CompoundSet":
        """Returns a copy of this set"""
        return CompoundSet(self.db, self.ids)

    def shuffled(self) -> "CompoundSet":
        """Returns a randomised copy of this set"""
        copy = self.copy()
        copy.shuffle()
        return copy

    def pop(self) -> Compound:
        """Pop the last compound in this set"""
        c_id = self.pop_id()
        return self.db.get_compound(id=c_id)

    def pop_id(self) -> int:
        """Pop the last compound id in this set"""
        return self._indices.pop()

    def shuffle(self) -> None:
        """Randomises the order of compounds in this set"""
        from random import shuffle

        shuffle(self._indices)

    def get_df(
        self,
        mol: bool = False,
        metadata: bool = False,
        poses: bool = False,
        count_by_target: bool = False,
        num_reactant: bool = False,
        num_reactions: bool = False,
        bases: bool = False,
        elabs: bool = False,
        tags: bool = False,
        num_poses: bool = False,
        num_routes: bool = False,
        alias: bool = True,
        **kwargs,
    ) -> "DataFrame":
        """Get a DataFrame representation of this set

        :param mol: include ``rdkit.Chem.Mol`` in output (Default value = False)
        :param metadata: include metadata in output (Default value = False)
        :param poses: include poses in output (Default value = False)
        :param count_by_target: count poses by target (Default value = False)
        :param num_reactant: include num_reactant column
        :param num_reactions: include num_reactions column
        :param bases: include bases column
        :param elabs: include elabs column
        :param tags: include tags column

        """

        from pandas import DataFrame

        data = []

        for comp in mrich.track(self, prefix="Creating compound dataframe"):
            d = comp.get_dict(
                mol=mol,
                metadata=metadata,
                count_by_target=count_by_target,
                poses=poses,
                num_reactant=num_reactant,
                num_reactions=num_reactions,
                bases=bases,
                elabs=elabs,
                tags=tags,
                **kwargs,
            )
            data.append(d)

        if num_poses:

            lookup = self.id_num_poses_dict

            for d in data:
                comp_id = d["id"]
                d["num_poses"] = lookup.get(comp_id, 0)

        df = DataFrame(data)

        if not alias:
            df.drop(columns=["alias"], inplace=True)

        return df

    def get_quoted(
        self,
        *,
        supplier: str = "any",
    ) -> "CompoundSet":
        """Get all member compounds that have a quote from given supplier

        :param supplier: supplier name (Default value = 'any')

        """

        if supplier == "any":
            key = f"quote_compound IN {self.str_ids}"
        else:
            key = f'quote_compound IN {self.str_ids} AND quote_supplier = "{supplier}"'

        ids = self.db.select_where(
            table="quote",
            query="DISTINCT quote_compound",
            key=key,
            multiple=True,
        )

        ids = [i for i, in ids]
        return CompoundSet(self.db, ids)

    def get_unquoted(
        self,
        *,
        supplier: str = "any",
    ) -> "CompoundSet":
        """Get all member compounds that do not have a quote from given supplier

        :param supplier: supplier name (Default value = 'any')

        """

        quoted = self.get_quoted(supplier=supplier)
        return self - quoted

    def get_dict(self) -> dict:
        """Get a dictionary object with all serialisable data needed to reconstruct this set"""
        return dict(db=str(self.db.path.resolve()), indices=self.indices)

    def write_smiles_csv(self, file: str) -> None:
        """Write a CSV of the smiles contained in this set to a file

        :param file: path of the CSV file

        """
        from pandas import DataFrame

        records = self.db.select_where(
            table=self.table,
            query="compound_id, compound_smiles",
            key=f"compound_id IN {self.str_ids}",
            multiple=True,
        )

        data = [dict(id=id, smiles=smiles) for id, smiles in records]
        df = DataFrame(data)
        mrich.writing(file)
        df.to_csv(file, index=False)

    def write_postera_csv(
        self,
        file,
        *,
        supplier: str = "Enamine",
        prefix: str = "fragment",
    ) -> None:
        """Write a CSV formatted for upload to Postera's Manifold

        :param file: path of the CSV file
        :param supplier: supplier to use for quotes, (Default value = 'Enamine')
        :param prefix: prefix to metadata columns, (Default value = 'fragment')

        """

        from datetime import date as dt
        from pandas import DataFrame

        if prefix:
            prefix = f"{prefix}_"

        data = []

        for c in mrich.track(self, prefix="Creating DataFrame"):

            # get props
            smiles = c.smiles
            tags = c.tags
            metadata = c.metadata
            poses = c.poses
            base = c.base

            # method
            assert len(tags) == 1, c
            method = tags[0]

            # date
            date = dt.today()

            # author
            assert "author" in metadata, c
            author = metadata["author"]

            match len(poses):
                case 1:
                    pose = poses[0]
                case 0:
                    mrich.warning(f"{c} has no poses")
                    assert base
                    pose = base.poses[0]
                case _:
                    mrich.warning(f"{c} has multiple poses")
                    pose = poses[0]

            # extract inspirations
            inspirations = pose.inspirations
            inspiration_names = ",".join(inspirations.names)
            inspiration_smiles = ".".join(inspirations.smiles)

            # quote info
            quotes = c.get_quotes(supplier=supplier)
            assert len(quotes) == 1, c
            quote = quotes[0]
            catalog_id = quote.entry
            catalog_price = quote.price
            catalog_lead_time = quote.lead_time

            # hippo string
            hippo_str = f"compound={c.id}, pose={pose.id}"

            # create row
            data.append(
                {
                    "SMILES": smiles,
                    f"{prefix}HIPPO_IDs": hippo_str,
                    f"{prefix}method": method,
                    f"{prefix}export_date": date,
                    f"{prefix}author": author,
                    f"{prefix}inspiration_names": inspiration_names,
                    f"{prefix}inspiration_SMILES": inspiration_smiles,
                    f"{prefix}supplier": supplier,
                    f"{prefix}supplier_catalogue": quote.catalogue,
                    f"{prefix}supplier_ID": catalog_id,
                    f"{prefix}supplier_price": catalog_price,
                    f"{prefix}supplier_lead_time": catalog_lead_time,
                }
            )

        df = DataFrame(data)

        mrich.writing(file)
        df.to_csv(file, index=False)

        return df

    def write_CAR_csv(
        self,
        file: "str | Path",
        amount: float = 1,  # in mg
        return_df: bool = False,
        # pick_cheapest: bool = False,
        quoted_only: bool = False,
        get_ingredient_quotes: bool = True,
        **kwargs,
    ) -> "DataFrame | None":
        """List of reactions for CAR

        Columns:

        * target-name
        * no-steps
        * concentration = None
        * amount-required
        * batch-tag

        per reaction

        * reactant-1-1
        * reactant-2-1
        * reaction-product-smiles-1
        * reaction-name-1
        * reaction-recipe-1
        * reaction-groupby-column-1

        :param file: output file
        :param amount: amount of each product in `mg`
        :param quoted_only: only choose reactants that have quotes
        :param supplier: only choose reactants that have quotes from this supplier
        :param kwargs: passed to :meth:`.Recipe.from_reaction`
        :param return_df: return a `DataFrame` (Default value = False)

        """

        from pathlib import Path
        from pandas import DataFrame
        from .recipe import Recipe

        file = str(Path(file).resolve())

        rows = []

        for r_id in mrich.track(self.reaction_ids, prefix="Solving compound recipes"):

            reaction = self.db.get_reaction(id=r_id)

            recipes = Recipe.from_reaction(
                reaction,
                amount=amount,
                pick_cheapest=False,
                quoted_only=quoted_only,
                get_ingredient_quotes=get_ingredient_quotes,
                **kwargs,
            )

            for sub_recipe in recipes:

                product = sub_recipe.product

                row = {
                    "target-names": str(product.compound),
                    "no-steps": 0,
                    "concentration-required-mM": None,
                    "amount-required-uL": None,
                    "batch-tag": None,
                }

                for i, reaction in enumerate(sub_recipe.reactions):

                    i = i + 1

                    row["no-steps"] += 1

                    match len(reaction.reactants):
                        case 1:
                            row[f"reactant-1-{i}"] = reaction.reactants[0].smiles
                            row[f"reactant-2-{i}"] = None
                        case 2:
                            row[f"reactant-1-{i}"] = reaction.reactants[0].smiles
                            row[f"reactant-2-{i}"] = reaction.reactants[1].smiles
                        case _:
                            raise NotImplementedError(
                                f"Unsupported number of reactants for {reaction=}: {len(reaction.reactants)}"
                            )

                    row[f"reaction-product-smiles-{i}"] = reaction.product.smiles
                    row[f"reaction-name-{i}"] = reaction.type
                    row[f"reaction-recipe-{i}"] = None
                    row[f"reaction-groupby-column-{i}"] = None
                    # row[f'reaction-id-{i}'] = int(reaction.id)

                rows.append(row)

        df = DataFrame(rows)

        df = df.convert_dtypes()

        for n_steps in set(df["no-steps"]):
            subset = df[df["no-steps"] == n_steps]
            this_file = file.replace(".csv", f"_{n_steps}steps.csv")
            mrich.writing(this_file)
            subset.to_csv(this_file, index=False)

        mrich.writing(file)
        df.to_csv(file, index=False)

        if return_df:
            return df

    def add_tag(
        self,
        tag: str,
    ) -> None:
        """Add this tag to every member of the set"""

        assert isinstance(tag, str)

        for i in self.indices:
            self.db.insert_tag(name=tag, compound=i, commit=False)

        mrich.print(f'Tagged {self} w/ "{tag}"')

        self.db.commit()

    ### DUNDERS

    def __len__(self) -> int:
        """The number of compounds in this set"""
        return len(self.indices)

    def __iter__(self):
        """Iterate through compounds in this set"""
        return iter(self.db.get_compound(id=i) for i in self.indices)

    def __getitem__(
        self,
        key: int | slice,
    ) -> "Compound | CompoundSet":
        """Get compounds or subsets thereof from this set

        :param key: integer index or slice of indices

        """
        match key:
            case int():
                index = self.indices[key]
                return self.db.get_compound(id=index)

            case slice():
                indices = self.indices[key]
                return CompoundSet(self.db, indices)

            case _:
                raise NotImplementedError

    def __sub__(
        self,
        other: "CompoundSet | IngredientSet",
    ) -> "CompoundSet":
        """Subtract a :class:`.Compound` object or ID to this set, or subtract multiple at once when ``other`` is a :class:`.CompoundSet` or :class:`.IngredientSet`"""

        match other:

            case CompoundSet():
                ids = set(self.ids) - set(other.ids)
                return CompoundSet(self.db, ids)

            case IngredientSet():
                mrich.warning(
                    "Subtracting IngredientSet from CompoundSet. Ignoring quote/amount data"
                )
                ids = set(self.ids) - set([int(i) for i in other.compound_ids])
                return CompoundSet(self.db, ids)

            case _:
                raise NotImplementedError

    def __add__(
        self,
        other: "Compound | CompoundSet | IngredientSet | int",
    ) -> "CompoundSet":
        """Add a :class:`.Compound` object or ID to this set, or add multiple at once when ``other`` is a :class:`.CompoundSet` or :class:`.IngredientSet`"""

        match other:

            case Compound():
                return self.add(other)

            case int():
                return self.add(other)

            case CompoundSet():
                ids = set(self.ids) | set(other.ids)
                return CompoundSet(self.db, ids)

            case IngredientSet():
                ids = set(self.ids) | set(other.compound_ids)
                return CompoundSet(self.db, ids)

            case _:
                raise NotImplementedError

    def __xor__(self, other: "CompoundSet"):
        """Exclusive OR set operation, returns all compounds in either set but not both"""

        match other:

            case CompoundSet():
                ids = set(self.ids) ^ set(other.ids)
                return CompoundSet(self.db, ids)

            case _:
                raise NotImplementedError

    def __str__(self) -> str:
        """Unformatted string representation"""

        if self.name:
            s = f"{self.name}: "
        else:
            s = ""

        s += "{" f"C × {len(self)}" "}"

        return s

    def __repr__(self) -> str:
        """ANSI ormatted string representation"""
        return f"{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}"

    def __rich__(self) -> str:
        """Representation for mrich"""
        return f"[bold underline]{self}"

    def __contains__(self, other: Compound | Ingredient | int):
        """Check if compound or ingredient is a member of this set"""
        match other:
            case Compound():
                id = other.id
            case Ingredient():
                id = other.compound_id
            case int():
                id = other

        return id in set(self.ids)


class IngredientSet:
    """An :class:`.Ingredient` is a :class:`.Compound` with a fixed quanitity and an attached quote, the :class:`.IngredientSet` is a object representing multiple ingredients.

    .. attention::

            :class:`.IngredientSet` objects should not be created directly. Instead they are returned by several methods when working with :doc:`quoting` and :doc:`rgen`.

    Selecting ingredients in the set
    ================================

    The :class:`.IngredientSet` can be indexed like a Python list:

    ::

            ingredient = ingredient_set[0] # first ingredient

    To get the ingredient for a specific :class:`.Compound` ID:

    ::

            ingredient = ingredient_set(compound_id=13)

    """

    _columns = [
        "compound_id",
        "amount",
        "quote_id",
        "supplier",
        "max_lead_time",
        "quoted_amount",
    ]

    def __init__(
        self,
        db: "Database",
        ingredients: "None | list[Ingredient]" = None,
        supplier: str | list | None = None,
        debug: bool = False,
    ) -> None:

        from pandas import DataFrame

        ingredients = ingredients or []

        self._db = db

        self._data = DataFrame(columns=self._columns, dtype=object)

        if debug:
            mrich.debug(self._data)

        self._supplier = supplier

        for ingredient in ingredients:
            self.add(ingredient)

        for col in self._columns:
            assert col in self._data.columns, f"{col} not in df.columns"

        if debug:
            mrich.debug(self._data)

    @classmethod
    def from_ingredient_df(
        cls,
        db: "Database",
        df: "DataFrame",
        supplier: str | list | None = None,
    ) -> "IngredientSet":
        """Create an :class:`.IngredientSet` from a DataFrame

        :param db: HIPPO Database
        :param df: DataFrame of Ingredients
        :param supplier: supplier to use for all quoting, (Default value = None)

        """
        # from numpy import nan
        self = cls.__new__(cls)

        for col in cls._columns:
            if col not in df.columns:
                raise Exception(f"{col} not in df.columns")
                df[col] = None

        self._db = db
        self._data = df.copy()
        self._supplier = supplier

        return self

    @classmethod
    def from_json(
        cls,
        db: "Database",
        path: None | str,
        supplier: str | list | None = None,
        data: None | dict = None,
    ) -> "IngredientSet":
        """Create an :class:`.IngredientSet` from JSON data or a JSON file

        :param db: HIPPO Database
        :param path: path to JSON data (can be ``None`` if ``data`` provided)
        :param supplier: supplier to use for all quoting, (Default value = ``None``)
        :param data: optional JSON data to parse, (Default value = ``None``)

        """

        if not data:
            import json

            data = json.load(open(path, "rt"))

        from pandas import DataFrame

        df = DataFrame(columns=cls._columns, dtype=object)

        for col in cls._columns:
            df[col] = data[col]

        return cls.from_ingredient_df(db=db, df=df, supplier=supplier)

    @classmethod
    def from_ingredient_dicts(
        cls,
        db: "Database",
        dicts: list[dict],
        supplier: str | list | None = None,
    ) -> "IngredientSet":
        """Create an :class:`.IngredientSet` from :class:`.Ingredient` dictionaries

        :param db: HIPPO Database
        :param dicts: List of individual ingredient dictionaries
        :param supplier: supplier to use for all quoting, (Default value = ``None``)

        """
        from pandas import DataFrame

        df = DataFrame(dicts, dtype=object)
        return cls.from_ingredient_df(db=db, df=df, supplier=supplier)

    @classmethod
    def from_compounds(
        cls,
        compounds: "CompoundSet | None" = None,
        ids: list[int] | None = None,
        db: "Database | None" = None,
        amount: float | list[float] = 1,
        supplier: str | list | None = None,
    ) -> "IngredientSet":
        """Create an :class:`.IngredientSet` from a :class:`.CompoundSet` or IDs

        :param compounds: :class:`.CompoundSet` to use, if ``None`` must provide ``ids`` and ``db`` (Default value = None)
        :param ids: Compound IDs (Default value = None)
        :param db: HIPPO Database (Default value = None)
        :param amount: Amount(s) in ``mg`` (Default value = 1)
        :param supplier: supplier to use for all quoting, (Default value = ``None``)

        """

        from pandas import DataFrame

        if not ids:
            ids = compounds.ids

        if not db:
            db = compounds.db

        df = DataFrame(
            dict(
                compound_id=ids,
                amount=amount,
                quote_id=None,
                supplier=supplier,
                max_lead_time=None,
                quoted_amount=None,
            ),
            dtype=object,
        )

        return cls.from_ingredient_df(db, df)

    ### PROPERTIES

    @property
    def df(self) -> "DataFrame":
        """Access the raw DataFrame"""
        return self._data

    @property
    def db(self) -> "Database":
        """Linked HIPPO Database"""
        return self._db

    @property
    def price_df(self) -> "DataFrame":
        """DataFrame including prices"""
        df = self.df.copy()
        df["price"] = [i.price for i in self]
        return df

    @property
    def price(self) -> "Price":
        """Total price of these ingredients"""
        return self.get_price()

    @property
    def supplier(self) -> str | list[str]:
        """Supplier(s)"""
        return self._supplier

    @supplier.setter
    def supplier(self, s):

        if isinstance(s, list) or isinstance(s, tuple):
            for x in s:
                assert isinstance(x, str)
        else:
            assert isinstance(s, str)

        self._supplier = s
        self.df["supplier"] = [s] * len(self)

    @property
    def smiles(self) -> list[str]:
        """SMILES for all ingredients"""
        compound_ids = list(self.df["compound_id"])
        result = self.db.select_where(
            query="compound_smiles",
            table="compound",
            key=f"compound_id in {tuple(compound_ids)}",
            multiple=True,
        )
        return [q for q, in result]

    @property
    def inchikeys(self) -> list[str]:
        """InChI-keys for all ingredients"""
        compound_ids = list(self.df["compound_id"])
        result = self.db.select_where(
            query="compound_inchikey",
            table="compound",
            key=f"compound_id in {tuple(compound_ids)}",
            multiple=True,
        )
        return [q for q, in result]

    @property
    def compound_ids(self) -> list[int]:
        """Compound IDs for all ingredients"""
        return list(self.df["compound_id"].values)

    @property
    def ids(self) -> list[int]:
        """Compound IDs for all ingredients"""
        return self.compound_ids

    @property
    def id_amount_pairs(self) -> list[tuple]:
        """Get a list of compound ID and amount pairs"""
        return [
            (id, amount) for id, amount in self.df[["compound_id", "amount"]].values
        ]

    @property
    def str_compound_ids(self) -> str:
        """Return an SQL formatted tuple string of the :class:`.Compound` IDs"""
        return str(tuple(self.df["compound_id"].values)).replace(",)", ")")

    @property
    def compounds(self) -> "CompoundSet":
        """:class:`.CompoundSet` of all compounds in this set"""
        return CompoundSet(self.db, self.compound_ids)

    ### METHODS

    def get_price(
        self, supplier: str | list[str] = None, none: str = "error"
    ) -> "Price":
        """Calculate the price with a given supplier

        :param supplier: supplier to use for all quoting, (Default value = ``None``)

        """

        from .price import Price

        pairs = {i: q for i, q in enumerate(self.df["quote_id"])}

        quote_ids = [q for q in pairs.values() if q is not None and not isnan(q)]

        if quote_ids:

            quote_id_str = str(tuple(quote_ids)).replace(",)", ")")

            if supplier:
                result = self.db.select_where(
                    query="quote_price, quote_currency",
                    table="quote",
                    key=f'quote_id in {quote_id_str} AND quote_supplier = "{supplier}"',
                    multiple=True,
                    none=none,
                )
            else:
                result = self.db.select_where(
                    query="quote_price, quote_currency",
                    table="quote",
                    key=f"quote_id in {quote_id_str}",
                    multiple=True,
                    none=none,
                )

            prices = [Price(a, b) for a, b in result]
            quoted = sum(prices, Price.null())

        else:

            quoted = Price.null()

        unquoted = [i for i, q in pairs.items() if q is None or isnan(q)]

        unquoted_price = Price.null()

        for i in unquoted:

            ingredient = self[i]

            p = ingredient.price

            unquoted_price += p

            quote = ingredient.quote

            assert quote, f"NULL Quote: {ingredient=}"

            self.df.loc[i, "quote_id"] = quote.id

            assert quote.amount

            self.df.loc[i, "quoted_amount"] = quote.amount

        return quoted + unquoted_price

    def interactive(self, **kwargs) -> None:
        """Wrapper for :meth:`.CompoundSet.interactive`"""
        self.compounds.interactive(**kwargs)

    def add(
        self,
        ingredient: "Ingredient | None" = None,
        *,
        compound_id: int | None = None,
        amount: float | None = None,
        quote_id: int | None = None,
        supplier: str | list[str] | None = None,
        max_lead_time: float | None = None,
        quoted_amount: float | None = None,
        debug: bool = False,
    ) -> None:
        """Add an

        :param ingredient: :class:`.Ingredient` to be added, if ``None`` must specify other parameters, (Default value = None)
        :param compound_id: :class:`.Compound` ID (Default value = None)
        :param amount: amount in ``mg`` (Default value = None)
        :param quote_id: :class:`.Quote` ID (Default value = None)
        :param supplier: supplier name string or list (Default value = None)
        :param max_lead_time: maximum lead-time for quoting (in days) (Default value = None)
        :param quoted_amount: amount of associated :class:`.Quote` (Default value = None)
        :param debug: increase verbosity for debugging (Default value = False)

        """

        from pandas import DataFrame, concat

        if ingredient:
            assert ingredient._table == "ingredient"
            compound_id = ingredient.compound_id
            amount = ingredient.amount

            if (q := ingredient.quote) and not ingredient.quote_id:
                mrich.warning(f"Losing quote! {ingredient.quote=}")

            supplier = ingredient.supplier
            max_lead_time = ingredient.max_lead_time

            if q is None:
                quote_id = None
                quoted_amount = None
            else:
                quote_id = q.id
                quoted_amount = q.amount

        else:
            assert compound_id
            assert amount

        if quote_id:
            # if not quoted_amount:
            #     mrich.warning(f'Requoting C{compound_id}...')

            assert quoted_amount

        supplier = self.supplier

        if self._data.empty:
            addition = DataFrame(
                [
                    dict(
                        compound_id=compound_id,
                        amount=amount,
                        quote_id=quote_id,
                        supplier=supplier,
                        max_lead_time=max_lead_time,
                        quoted_amount=quoted_amount,
                    )
                ],
                dtype=object,
            )
            self._data = addition

        else:

            if compound_id in self._data["compound_id"].values:
                index = self._data.index[
                    self._data["compound_id"] == compound_id
                ].tolist()[0]
                self._data.loc[index, "amount"] += amount

                # discard if the quote is no longer valid
                if (a := self.df.loc[index, "quoted_amount"]) and a < self.df.loc[
                    index, "amount"
                ]:
                    self._data.loc[index, "quote_id"] = None
                    self._data.loc[index, "quoted_amount"] = None

                if debug and supplier:
                    mrich.debug("Adding to existing ingredient")
                    mrich.debug(f'{self._data.loc[index, "supplier"]=}')
                    mrich.debug(f"{supplier=}")

            else:
                # from numpy import nan
                addition = DataFrame(
                    [
                        dict(
                            compound_id=compound_id,
                            amount=amount,
                            quote_id=quote_id,
                            supplier=supplier,
                            max_lead_time=max_lead_time,
                            quoted_amount=quoted_amount,
                        )
                    ],
                    dtype=object,
                )

                self._data = concat(
                    [self._data, addition], ignore_index=True, join="inner"
                )

                if debug:
                    mrich.out(addition)

    def _get_ingredient(
        self,
        series,
    ) -> "Ingredient":
        """Get ingredient from one of the DataFrame rows"""

        q_id = series["quote_id"]

        if isinstance(q_id, float) and isnan(q_id):
            q_id = None

        return Ingredient(
            db=self._db,
            compound=series["compound_id"],
            amount=series["amount"],
            quote=q_id,
            supplier=series["supplier"],
            max_lead_time=series["max_lead_time"],
        )

    def copy(self) -> "IngredientSet":
        """Return a copy of this :class:`.IngredientSet`"""
        return IngredientSet.from_ingredient_df(
            self.db, self.df, supplier=self.supplier
        )

    def draw(self) -> None:
        """Wrapper for :meth:`.CompoundSet.draw`"""
        self.compounds.draw()

    def set_amounts(
        self,
        amount: float | list[float],
    ) -> None:
        """Set the amount(s) for all ingredients in this set, and update quotes

        :param amount: amount in ``mg``

        """

        self.df["amount"] = amount

        # if amounts are modified the quotes should be cleared
        self.df["quote_id"] = None

        assert all(self.df["supplier"].isna()) and all(self.df["max_lead_time"].isna())

        # update quotes
        pairs = self.db.execute(
            f"""
            WITH matching_quotes AS (
                SELECT quote_id, quote_compound, MIN(quote_price) FROM quote
                WHERE quote_compound IN {self.str_compound_ids}
                AND quote_amount >= {amount}
                GROUP BY quote_compound
            )
            SELECT compound_id, quote_id FROM compound
            LEFT JOIN matching_quotes ON quote_compound = compound_id
            WHERE compound_id IN {self.str_compound_ids}
        """
        ).fetchall()

        for compound_id, quote_id in pairs:
            match = self.df.index[self.df["compound_id"] == compound_id][0]
            self.df.loc[match, "quote_id"] = quote_id

    def get_dict(self, data_orient: str = "list") -> dict:
        """Get serialisable dictionary

        :param data_orient: passed to ``pandas.DataFrame.to_dict`` (Default value = 'list')

        """
        return dict(
            db=str(self.db),
            supplier=self.supplier,
            data=self.df.to_dict(orient=data_orient),
        )

    ### DUNDERS

    def __len__(self):
        """The number of ingredients in this set"""
        return len(self._data)

    def __str__(self) -> str:
        """Unformatted string representation"""
        return "{" f"Ingredient × {len(self)}" "}"

    def __repr__(self) -> str:
        """ANSI ormatted string representation"""
        return f"{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}"

    def __rich__(self) -> str:
        """Representation for mrich"""
        return f"[bold underline]{self}"

    def __add__(self, other):
        """Add another  :class:`.IngredientSet` this set"""

        for i, row in other._data.iterrows():
            self.add(
                compound_id=row.compound_id,
                amount=row.amount,
                quote_id=row.quote_id,
                supplier=row.supplier,
                max_lead_time=row.max_lead_time,
                quoted_amount=row.quoted_amount,
            )

        return self

    def __getitem__(self, key):
        match key:
            case int():
                series = self.df.loc[key]
                return self._get_ingredient(series)

            case _:
                raise NotImplementedError

    def __iter__(self):
        return iter(self._get_ingredient(s) for i, s in self.df.iterrows())

    def __call__(
        self,
        *,
        compound_id: int | None = None,
        tag: str | None = None,
    ) -> "IngredientSet | Ingredient | CompoundSet":

        if compound_id:

            # get the ingredient with the matching compound ID
            matches = self.df[self.df["compound_id"] == compound_id]

            if len(matches) == 0:
                return None

            elif len(matches) != 1:

                mrich.warning(f"Multiple ingredients in set with {compound_id=}")
                # print(matches)

                return IngredientSet(
                    self.db, [self._get_ingredient(s) for i, s in matches.iterrows()]
                )

            return self._get_ingredient(matches.iloc[0])

        # elif tag:
        #     return self.compounds(tag=tag)

        else:
            raise NotImplementedError

    def __getattr__(self, key: str):
        """For missing attributes try getting from associated :class:`.CompoundSet`"""
        return getattr(self.compounds, key)

    def __contains__(self, other: Compound | Ingredient | int):
        """Check if compound or ingredient is a member of this set"""
        match other:
            case Compound():
                id = other.id
            case Ingredient():
                id = other.compound_id
            case int():
                id = other

        return id in set(self.compound_ids)
