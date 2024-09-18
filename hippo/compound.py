import mcol

from rdkit import Chem

from .pose import Pose

# from .pset import PoseSet
from .tags import TagSet

# from .rset import ReactionSet
from .target import Target
from .quote import Quote

import logging

logger = logging.getLogger("HIPPO")


class Compound:
    """A :class:`.Compound` represents a ligand/small molecule with stereochemistry removed and no atomic coordinates. I.e. it represents the chemical structure. It's name is always an InChiKey. If a compound is an elaboration it can have a :meth:`.Compound.base` property which is another :class:`.Compound`. :class:`.Compound` objects are target-agnostic and can be linked to any number of catalogue entries (:class:`.Quote`) or synthetic pathways (:class:`.Reaction`).

    .. attention::

            :class:`.Compound` objects should not be created directly. Instead use :meth:`.HIPPO.register_compound` or :meth:`.HIPPO.compounds`. See :doc:`getting_started` and :doc:`insert_elaborations`.

    """

    _table = "compound"

    def __init__(
        self,
        animal: "HIPPO",
        db: "Database",
        id: int,
        inchikey: str,
        alias: str,
        smiles: str,
        mol: Chem.Mol | bytes | None = None,
        metadata: dict | None = None,
    ):

        # from compound table
        self._id = id
        self._inchikey = inchikey
        self._alias = alias
        self._smiles = smiles
        self._animal = animal
        self._bases = None
        self._elabs = None
        self._alias = alias
        self._tags = None
        self._metadata = metadata

        # computed properties
        self._num_heavy_atoms = None
        self._num_rings = None
        self._formula = None
        self._molecular_weight = None
        self._total_changes = db.total_changes

        if isinstance(mol, bytes):
            mol = Chem.Mol(mol)

        self._mol = mol

        self._db = db

    ### FACTORIES

    ### PROPERTIES

    @property
    def id(self) -> int:
        """Returns the compound's database ID"""
        return self._id

    @property
    def inchikey(self) -> str:
        """Returns the compound's InChiKey"""
        return self._inchikey

    @property
    def name(self) -> str:
        """Returns the compound's InChiKey"""
        return self._inchikey

    @property
    def smiles(self) -> str:
        """Returns the compound's (flattened) smiles"""
        return self._smiles

    @property
    def alias(self) -> str:
        """Returns the compound's alias"""
        return self._alias

    @alias.setter
    def alias(self, alias: str) -> None:
        """Set the compound's alias"""
        self.set_alias(alias)

    @property
    def mol(self) -> Chem.Mol:
        """Returns the compound's RDKit Molecule"""
        if self._mol is None:
            (mol,) = self.db.select_where(
                query="mol_to_binary_mol(compound_mol)",
                table="compound",
                key="id",
                value=self.id,
                multiple=False,
            )
            self._mol = Chem.Mol(mol)
        return self._mol

    @property
    def num_heavy_atoms(self) -> int:
        """Get the number of heavy atoms"""
        if self._num_heavy_atoms is None:
            self._num_heavy_atoms = self.db.get_compound_computed_property(
                "num_heavy_atoms", self.id
            )
        return self._num_heavy_atoms

    @property
    def molecular_weight(self) -> float:
        """Get the molecular weight"""
        if self._molecular_weight is None:
            self._molecular_weight = self.db.get_compound_computed_property(
                "molecular_weight", self.id
            )
        return self._molecular_weight

    @property
    def num_rings(self) -> int:
        """Get the number of rings"""
        if self._num_rings is None:
            self._num_rings = self.db.get_compound_computed_property(
                "num_rings", self.id
            )
        return self._num_rings

    @property
    def formula(self) -> str:
        """Get the chemical formula"""
        if self._formula is None:
            self._formula = self.db.get_compound_computed_property("formula", self.id)
        return self._formula

    @property
    def atomtype_dict(self) -> dict[str, int]:
        """Get a dictionary with atomtypes as keys and corresponding quantities/counts as values."""
        from molparse.atomtypes import formula_to_atomtype_dict

        return formula_to_atomtype_dict(self.formula)

    @property
    def num_atoms_added(self) -> int:
        """Calculate the number of atoms added relative to the base compound"""
        assert (b_id := self._base_id), f"{self} has no base defined"
        n_e = self.num_heavy_atoms
        n_b = self.db.get_compound_computed_property("num_heavy_atoms", b_id)
        return n_e - n_b

    @property
    def metadata(self) -> "MetaData":
        """Returns the compound's metadata dict"""
        if self._metadata is None:
            self._metadata = self.db.get_metadata(table="compound", id=self.id)
        return self._metadata

    @property
    def db(self) -> "Database":
        """Returns a pointer to the parent database"""
        return self._db

    @property
    def tags(self) -> TagSet:
        """Returns the compound's tags"""
        if not self._tags:
            self._tags = self.get_tags()
        return self._tags

    @property
    def poses(self) -> "PoseSet":
        """Returns the compound's poses"""
        return self.get_poses()

    @property
    def best_placed_pose(self) -> Pose:
        """Returns the compound's pose with the lowest distance score"""
        return self.poses.best_placed_pose

    @property
    def num_poses(self) -> int:
        """Returns the number of associated poses"""
        return self.db.count_where(table="pose", key="compound", value=self.id)

    @property
    def num_reactions(self) -> int:
        """Returns the number of associated reactions (product)"""
        return self.db.count_where(table="reaction", key="product", value=self.id)

    @property
    def num_reactant(self) -> int:
        """Returns the number of associated reactions (reactant)"""
        return self.db.count_where(table="reactant", key="compound", value=self.id)

    @property
    def bases(self):
        """Returns the base compound for this elaboration"""
        if self._bases is None or self._db_changed:
            ids = self.get_base_ids()
            if not ids:
                self._bases = None
            else:
                from .cset import CompoundSet

                self._bases = CompoundSet(
                    self.db, ids, name=f"base scaffolds of {self}"
                )
                self._total_changes = self.db.total_changes
        return self._bases

    @property
    def elabs(self):
        """Returns the base compound for this elaboration"""
        if self._elabs is None or self._db_changed:
            ids = self.get_superstructure_ids()
            if not ids:
                self._elabs = None
            else:
                from .cset import CompoundSet

                self._elabs = CompoundSet(self.db, ids, name=f"elaborations of {self}")
                self._total_changes = self.db.total_changes
        return self._elabs

    @property
    def reactions(self) -> "ReactionSet":
        """Returns the reactions resulting in this compound"""
        return self.get_reactions(none=False)

    @property
    def reaction(self) -> "Reaction":
        """Returns the reaction resulting in this compound (will return first if multiple, with a warning)"""
        reactions = self.reactions
        match len(reactions):
            case 0:
                logger.warning(f"{self} has no reactions")
                return None
            case 1:
                logger.warning(f"{self} has multiple reactions, returning first")
            case _:
                pass

        return reactions[0]

    @property
    def dict(self) -> dict:
        """Returns a dictionary of this compound. See :meth:`.Compound.get_dict`"""
        return self.get_dict()

    @property
    def is_base(self) -> bool:
        """Is this Compound the basis for any elaborations?"""
        return bool(
            self.db.select_where(
                query="1",
                table="scaffold",
                key="base",
                value=self.id,
                multiple=False,
                none="quiet",
            )
        )

    @property
    def is_elab(self) -> bool:
        """Is this Compound the based on any other compound?"""
        return bool(
            self.db.select_where(
                query="1",
                table="scaffold",
                key="superstructure",
                value=self.id,
                multiple=False,
                none="quiet",
            )
        )

    @property
    def is_product(self) -> bool:
        """Is this Compound a product of at least one reaction"""
        return bool(self.get_reactions(none=False))

    @property
    def table(self):
        """Returns the name of the :class:`.Database` table"""
        return self._table

    @property
    def _db_changed(self) -> bool:
        """Has the database changed?"""
        return self._total_changes != self.db.total_changes

    ### METHODS

    def add_stock(
        self,
        amount: float,
        *,
        purity: float | None = None,
        entry: str | None = None,
        location: str | None = None,
        return_quote: bool = True,
    ) -> int | Quote:
        """Register a certain quantity of compound stock in the Database.

        :param amount: Amount in ``mg``
        :param purity: Purity fraction ``0 < purity <= 1``, defaults to ``None``
        :param location: String describing where this stock is located, defaults to ``None``
        :param return_quote: If ``True`` a :class:`.Quote` object is returned instead of its ID, defaults to ``True``
        :returns: The inserted :class:`.Quote` object or ID (see ``return_quote``)
        """

        assert amount

        # search for existing in stock quotes
        existing = self.get_quotes(supplier="Stock", df=False)

        # supersede old in stock records
        if existing:
            delete = set()
            not_deleted = 0
            for quote in existing:

                if any(
                    [
                        quote.entry != entry,
                        quote.purity != purity,
                        quote.catalogue != location,
                    ]
                ):
                    not_deleted += 1
                    continue

                delete.add(quote.id)

            delete_str = str(tuple(delete)).replace(",)", ")")

            self.db.delete_where(table="quote", key=f"quote_id IN {delete_str}")

            if delete:
                logger.warning(f"Removed {len(delete)} existing In-Stock Quotes")

            if not_deleted:
                logger.warning(
                    f"Did not remove {not_deleted} existing In-Stock Quotes with differing entry/purity/location"
                )

        # insert the new quote
        quote_id = self.db.insert_quote(
            compound=self.id,
            price=0,
            lead_time=0,
            currency=None,
            supplier="Stock",
            catalogue=location,
            entry=entry,
            amount=amount,
            purity=purity,
        )

        if return_quote:
            self.db.get_quote(id=quote_id)
        else:
            return quote_id

    def get_tags(self) -> "TagSet":
        """Get the tags assigned to this compound"""
        tags = self.db.select_where(
            query="tag_name",
            table="tag",
            key="compound",
            value=self.id,
            multiple=True,
            none="quiet",
        )
        return TagSet(self, {t[0] for t in tags}, commit=False)

    def get_quotes(
        self,
        min_amount: float | None = None,
        supplier: str | None = None,
        max_lead_time: float | None = None,
        none: str = "quiet",
        pick_cheapest: bool = False,
        df: bool = False,
    ) -> list["Quote"]:
        """Get all quotes associated to this compound

        :param min_amount: Only return quotes with amounts greater than this, defaults to ``None``
        :param supplier: Only return quotes with the given supplier, defaults to ``None``
        :param max_lead_time: Only return quotes with lead times less than this (in days), defaults to ``None``
        :param none: Define the behaviour when no quotes are found. Choose `error` to raise print an error.
        :param pick_cheapest: If ``True`` only the cheapest :class:`.Quote` is returned, defaults to ``False``
        :param df: Returns a ``DataFrame`` of the quoting data, defaults to ``False``
        :returns: List of :class:`.Quote` objects, ``DataFrame``, or single :class:`.Quote`. See ``pick_cheapest`` and ``df`` parameters

        """

        if not supplier:
            quote_ids = self.db.select_where(
                query="quote_id",
                table="quote",
                key="compound",
                value=self.id,
                multiple=True,
                none=none,
            )
        elif isinstance(supplier, str):
            quote_ids = self.db.select_where(
                query="quote_id",
                table="quote",
                key=f'quote_compound = {self.id} AND quote_supplier = "{supplier}"',
                multiple=True,
                none=none,
            )
        else:
            quote_ids = self.db.select_where(
                query="quote_id",
                table="quote",
                key=f'quote_compound = {self.id} AND quote_supplier IN {str(tuple(supplier)).replace(",)",")")}',
                multiple=True,
                none=none,
            )

        if quote_ids:
            quotes = [self.db.get_quote(id=q[0]) for q in quote_ids]
        else:
            return None

        if max_lead_time:
            quotes = [q for q in quotes if q.lead_time <= max_lead_time]

        if min_amount:
            suitable_quotes = [q for q in quotes if q.amount >= min_amount]

            if not suitable_quotes:
                logger.debug(f"No quote available with amount >= {min_amount} mg")
                quotes = [Quote.combination(min_amount, quotes)]

            else:
                quotes = suitable_quotes

        if pick_cheapest:
            return sorted(quotes, key=lambda x: x.price)[0]

        if df:
            from pandas import DataFrame

            return DataFrame([q.dict for q in quotes]).drop(columns="compound")

        return quotes

    def get_reactions(
        self,
        as_reactant: bool = False,
        permitted_reactions: "ReactionSet" = None,
        none: str = "error",
    ) -> "ReactionSet":
        """Get the associated :class:`.Reaction` objects. By default this function returns all reaction resulting in this :class:`.Compound` as a product, unless ``as_reactant`` is set to ``True``.

        :param as_reactant: Search for :class:`.Reaction` objects using this :class:`.Compound` as a reactant instead of a product, defaults to ``False``
        :param permitted_reactions: Provide a :class:`.ReactionSet` by which to filter the results
        :param none: Define the behaviour when no quotes are found. Choose `error` to raise print an error, defaults to ``'error'``
        """

        from .rset import ReactionSet

        if as_reactant:
            reaction_ids = self.db.select_where(
                query="reactant_reaction",
                table="reactant",
                key="compound",
                value=self.id,
                multiple=True,
                none=none,
            )
        else:
            reaction_ids = self.db.select_where(
                query="reaction_id",
                table="reaction",
                key="product",
                value=self.id,
                multiple=True,
                none=none,
            )

        reaction_ids = [q for q, in reaction_ids]

        if permitted_reactions:
            reaction_ids = [i for i in reaction_ids if i in permitted_reactions]

        return ReactionSet(self.db, reaction_ids)

    def get_poses(self) -> "PoseSet":
        """Get the associated :class:`.Pose` objects."""

        pose_ids = self.db.select_where(
            query="pose_id",
            table="pose",
            key="compound",
            value=self.id,
            multiple=True,
            none=False,
        )

        from .pset import PoseSet

        return PoseSet(self.db, [q[0] for q in pose_ids], name=f"{self}'s poses")

    def get_dict(
        self,
        *,
        mol: bool = True,
        metadata: bool = True,
        poses: bool = True,
        count_by_target: bool = False,
    ) -> "dict":
        """Returns a dictionary representing this :class:`.Compound`

        :param mol: Include a ``rdkit.Chem.Mol object``, defaults to ``True``
        :param metadata: Include metadata, defaults to ``True``
        :param poses: Include dictionaries of associated :class:`.Pose` objects, defaults to ``True``
        :param count_by_target: Include counts by protein :class:`.Target`, defaults to ``False``. Only applicable when ``count_by_target = True``.
        :returns: A dictionary
        """

        serialisable_fields = [
            "id",
            "alias",
            "inchikey",
            "smiles",
            "num_reactant",
            "num_reactions",
        ]

        data = {}
        for key in serialisable_fields:
            data[key] = getattr(self, key)

        if mol:
            try:
                data["mol"] = self.mol
            except InvalidMolError:
                data["mol"] = None

        if self.base:
            data["base"] = self.base.name
        else:
            data["base"] = None

        data["tags"] = self.tags

        if poses:

            poses = self.poses

            if poses:

                data["poses"] = poses.ids
                data["targets"] = poses.target_names

                if count_by_target:
                    target_ids = poses.target_ids

                    for target in self._animal.targets:
                        t_poses = poses(target=target.id) or []
                        data[f"#poses {target.name}"] = len(t_poses)

        if metadata and (metadict := self.metadata):
            for key in metadict:
                data[key] = metadict[key]

        return data

    def get_base_ids(self) -> list[int]:
        """Get a list of :class:`.Compound` ID's that this object is a superstructure of"""
        ids = self.db.select_where(
            table="scaffold",
            query="scaffold_base",
            key="superstructure",
            value=self.id,
            none="quiet",
            multiple=True,
        )
        if not ids:
            return None
        return [i for i, in ids]

    def get_superstructure_ids(self) -> list[int]:
        """Get a list of :class:`.Compound` ID's that this object is a substructure of"""
        ids = self.db.select_where(
            table="scaffold",
            query="scaffold_superstructure",
            key="base",
            value=self.id,
            none="quiet",
            multiple=True,
        )
        if not ids:
            return None
        return [i for i, in ids]

    def add_base(self, base: "Compound | int", commit: bool = True) -> None:
        """
        Add a base :class:`.Compound` this molecule is derived from.

        :param base: The base :class:`.Compound` or its ID.
        :param commit: Commit the changes to the :class:`.Database`, defaults to ``True``
        """

        if not isinstance(base, int):
            assert base._table == "compound"
            base = base.id
        self.db.insert_scaffold(base=base, superstructure=self.id)

    def set_alias(self, alias: str, commit=True) -> None:
        """
        Set this :class:`.Compound`'s alias.

        :param alias: The alias
        :param commit: Commit the changes to the :class:`.Database`, defaults to ``True``
        """

        assert isinstance(alias, str)
        self._alias = alias
        self.db.update(
            table="compound",
            id=self.id,
            key="compound_alias",
            value=alias,
            commit=commit,
        )

    def as_ingredient(
        self,
        amount: float,
        max_lead_time: float = None,
        supplier: str = None,
    ) -> "Ingredient":
        """Convert this compound into an :class:`.Ingredient` object with an associated amount (in ``mg``) and :class:`.Quote` if available.

        :param amount: Amount in ``mg``
        :param supplier: Only search for quotes with the given supplier, defaults to ``None``
        :param max_lead_time: Only search for quotes with lead times less than this (in days), defaults to ``None``
        """

        quote = self.get_quotes(
            pick_cheapest=True,
            min_amount=amount,
            max_lead_time=max_lead_time,
            supplier=supplier,
            none="quiet",
        )

        if not quote:
            quote = None

        return Ingredient(
            db=self.db,
            compound=self.id,
            amount=amount,
            quote=quote,
            supplier=supplier,
            max_lead_time=max_lead_time,
        )

    def draw(self, align_substructure: bool = False) -> None:
        """Display this compound (and its base if it has one)

        .. attention::

                This method is only intended for use within a Jupyter Notebook.

        :param align_substructure: Align the two drawing by their common substructure, defaults to ``False``
        """

        if base := self.base:

            from molparse.rdkit import draw_mcs

            data = {self.base.smiles: f"{base} (base)", self.smiles: str(self)}

            if len(data) == 2:

                drawing = draw_mcs(
                    data,
                    align_substructure=align_substructure,
                    show_mcs=False,
                    highlight_common=False,
                )
                display(drawing)

            else:
                logger.error(
                    f"Problem drawing {base.id=} vs {self.id=}, self referential?"
                )
                display(self.mol)

        else:
            display(self.mol)

    def classify(
        self,
        draw: bool = True,
    ) -> list[tuple[str, int]]:
        """
        Find RDKit Fragments within the compound molecule and draw them

        :param draw: Draw the annotated molecule, defaults to ``True``
        :returns: A list of tuples containing a descriptor (``str``) and count (``int``) pair
        """
        # from molparse.rdkit import classify_mol
        from molparse.rdkit.classify import classify_mol

        return classify_mol(self.mol, draw=draw)

    def murcko_scaffold(
        self,
        generic: bool = False,
    ):
        """Get the rdkit MurckoScaffold for this compound"""

        from rdkit.Chem.Scaffolds import MurckoScaffold

        scaffold = MurckoScaffold.GetScaffoldForMol(self.mol)

        if generic:
            scaffold = MurckoScaffold.MakeScaffoldGeneric(scaffold)

        return scaffold

    def summary(self, metadata: bool = True, draw: bool = True) -> None:
        """
        Print a summary of this compound

        :param metadata: Include metadata, defaults to ``True``
        :param draw: Include a 2D molecule drawing, defaults to ``True``
        """

        logger.header(repr(self))

        logger.var("inchikey", self.inchikey)
        logger.var("alias", self.alias)
        logger.var("smiles", self.smiles)
        logger.var("base", self.base)

        logger.var("is_base", self.is_base)
        logger.var("num_heavy_atoms", self.num_heavy_atoms)
        logger.var("num_rings", self.num_rings)
        logger.var("formula", self.formula)

        poses = self.poses
        logger.var("#poses", len(poses))
        if poses:
            logger.var("targets", poses.targets)

        logger.var("#reactions (product)", self.num_reactions)
        logger.var("#reactions (reactant)", self.num_reactant)

        logger.var("tags", self.tags)

        if metadata:
            logger.var("metadata", str(self.metadata))

        if draw:
            self.draw()

    def place(
        self,
        *,
        reference: Pose,
        inspirations: list[Pose] | None = None,
        max_ddG: float = 0.0,
        max_RMSD: float = 2.0,
        output_dir: str = "wictor_place",
        tags: list[str] = None,
        metadata: dict = None,
        overwrite: bool = False,
    ) -> Pose:
        """
        Generate a new pose for this compound using Fragmenstein.

        :param reference: Choose the :class:`.Pose` to use as the reference protein conformation
        :param inspirations: Choose the (virtual) hits to to define the ligand reference, defaults to the ``reference``'s inspirations
        :param max_ddG: Maximum ``ddG`` value permitted for a valid ligand conformation, defaults to ``0.0``
        :param max_RMSD: Maximum ``RMSD`` value permitted for a valid ligand conformation, defaults to ``2.0``
        :param output_dir: Output directory for Fragmenstein files, defaults to ``wictor_place``
        :param tags: Tags to assign to the created pose, defaults to ``[]``
        :param metadata: A dictionary of metadata to assign to this compound, defaults to ``{}``
        :param overwrite: Delete old poses, defaults to ``False``
        """

        from fragmenstein import Monster, Wictor
        from pathlib import Path

        tags = tags or []
        metadata = metadata or {}

        # get required data
        smiles = self.smiles

        inspirations = inspirations or reference.inspirations
        target = reference.target.name

        inspiration_mols = [c.mol for c in inspirations]
        protein_pdb_block = reference.protein_system.pdb_block_with_alt_sites

        # create the victor
        victor = Wictor(hits=inspiration_mols, pdb_block=protein_pdb_block)
        victor.work_path = output_dir
        victor.enable_stdout(logging.CRITICAL)

        # do the placement
        victor.place(smiles, long_name=self.name)

        # metadata
        metadata["ddG"] = (
            victor.energy_score["bound"]["total_score"]
            - victor.energy_score["unbound"]["total_score"]
        )
        metadata["RMSD"] = victor.mrmsd.mrmsd

        if metadata["ddG"] > max_ddG:
            return None

        if metadata["RMSD"] > max_RMSD:
            return None

        # register the pose
        pose = self._animal.register_pose(
            compound=self,
            target=target,
            path=Path(victor.work_path) / self.name / f"{self.name}.minimised.mol",
            inspirations=inspirations,
            reference=reference,
            tags=tags,
            metadata=metadata,
        )

        if overwrite:
            ids = [p.id for p in self.poses if p.id != pose.id]
            for i in ids:
                self.db.delete_where(table="pose", key="id", value=i)
            logger.success(f"Successfully posed {self} (and deleted old poses)")
        else:
            logger.success(f"Successfully posed {self}")

        return pose

    ### DUNDERS

    def __str__(self):
        return f"C{self.id}"

    def __repr__(self):
        return f'{mcol.bold}{mcol.underline}{self} "{self.name}"{mcol.unbold}{mcol.ununderline}'

    def __eq__(self, other):
        return self.id == other.id


class Ingredient:
    """An ingredient is a :class:`.Compound` with a fixed quanitity and an attached quote.

    .. image:: ../images/ingredient.png
              :width: 450
              :alt: Ingredient schema

    .. attention::

            :class:`.Ingredient` objects should not be created directly. Instead use :meth:`.Compound.as_ingredient`.
    """

    _table = "ingredient"

    def __init__(
        self, db, compound, amount, quote, max_lead_time=None, supplier=None
    ) -> "Ingredient":

        assert compound

        self._db = db

        # don't store inherited compound in memory until needed
        self._compound = None

        if isinstance(compound, Compound):
            self._compound_id = compound.id
            self._compound = None
        else:
            self._compound_id = compound

        if isinstance(quote, Quote):

            if id := quote.id:
                self._quote_id = quote.id
                self._quote = None

            else:
                self._quote_id = None
                self._quote = quote

        elif quote is None:
            self._quote_id = None
            self._quote = None

        else:
            self._quote_id = int(quote)
            self._quote = None

        self._amount = amount
        self._max_lead_time = max_lead_time
        self._supplier = supplier

    ### PROPERTIES

    @property
    def db(self) -> "Database":
        """Returns the parent :class:`.Database`"""
        return self._db

    @property
    def amount(self) -> float:
        """Returns the amount (in ``mg``)"""
        return self._amount

    @property
    def id(self) -> int:
        """Returns the ID of the associated :class:`.Compound`"""
        return self._compound_id

    @property
    def compound_id(self) -> int:
        """Returns the ID of the associated :class:`.Compound`"""
        return self._compound_id

    @property
    def quote_id(self) -> int:
        """Returns the ID of the associated :class:`.Quote`"""
        return self._quote_id

    @property
    def max_lead_time(self) -> float:
        """Returns the max_lead_time (in days) from the original quote query"""
        return self._max_lead_time

    @property
    def supplier(self) -> str:
        """Returns the supplier from the original quote query"""
        return self._supplier

    @amount.setter
    def amount(self, a) -> None:
        """Set the amount and fetch updated :class:`.Quote`s"""

        quote_id = self.get_cheapest_quote_id(
            min_amount=a,
            max_lead_time=self._max_lead_time,
            supplier=self._supplier,
            none="quiet",
        )

        self._quote_id = quote_id

        self._amount = a

    @property
    def compound(self) -> Compound:
        """Returns the associated :class:`.Compound`"""

        if not self._compound:
            self._compound = self.db.get_compound(id=self.compound_id)
        return self._compound

    @property
    def quote(self) -> Quote:
        """Returns the associated :class:`.Quote`"""

        if not self._quote:
            if q_id := self.quote_id:
                self._quote = self.db.get_quote(id=self.quote_id)

            else:
                q = self.compound.get_quotes(
                    pick_cheapest=True,
                    min_amount=self.amount,
                    max_lead_time=self.max_lead_time,
                    supplier=self.supplier,
                    none="error",
                )

                if not q:
                    return None

                self._quote = q
                self._quote_id = q.id

        return self._quote

    @property
    def compound_price_amount_str(self) -> str:
        """String representation including :class:`.Compound`, :class:`.Price`, and amount."""
        return f"{self} ({self.amount})"

    @property
    def smiles(self) -> str:
        """Returns the SMILES of the associated :class:`.Compound`"""
        return self.compound.smiles

    @property
    def price(self) -> "Price | None":
        """Returns the :class:`.Price` of the associated :class:`.Quote`"""
        if self.quote:
            return self.quote.price
        else:
            return None

    @property
    def lead_time(self) -> float | None:
        """Returns the lead time (in days) of the associated :class:`.Quote`"""
        if self.quote:
            return self.quote.lead_time
        else:
            return None

    ### METHODS

    def get_cheapest_quote_id(
        self,
        min_amount: float | None = None,
        supplier: str | None = None,
        max_lead_time: float | None = None,
        none: str = "quiet",
    ) -> int | None:
        """
        Query quotes associated to this ingredient, and return the cheapest

        :param min_amount: Only return quotes with amounts greater than this, defaults to ``None``
        :param supplier: Only return quotes with the given supplier, defaults to ``None``
        :param max_lead_time: Only return quotes with lead times less than this (in days), defaults to ``None``
        :param none: Define the behaviour when no quotes are found. Choose `error` to raise print an error.
        """

        supplier_str = f' AND quote_supplier IS "{supplier}"' if supplier else ""
        lead_time_str = (
            f" AND quote_lead_time <= {max_lead_time}" if max_lead_time else ""
        )
        key_str = f"quote_compound IS {self.compound_id} AND quote_amount >= {min_amount}{supplier_str}{lead_time_str} ORDER BY quote_price"

        result = self.db.select_where(
            query="quote_id", table="quote", key=key_str, multiple=False, none=none
        )

        if result:
            (quote_id,) = result
            return quote_id

        else:
            return None

    def get_quotes(self, **kwargs) -> list["Quote"]:
        """Wrapper for :meth:`.Compound.get_quotes()`"""
        return self.compound.get_quotes(**kwargs)

    ### DUNDERS

    def __str__(self) -> str:
        """Plain string representation"""
        return f"{self.amount:.2f}mg of C{self._compound_id}"

    def __repr__(self) -> str:
        """Formatted string representation"""
        return f"{mcol.bold}{mcol.underline}{str(self)}{mcol.unbold}{mcol.ununderline}"

    def __eq__(self, other) -> bool:
        """Equality operator"""

        if self.compound_id != other.compound_id:
            return False

        return self.amount == other.amount
