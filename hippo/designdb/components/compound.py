"""Compound and Ingredient components."""

import logging
from pathlib import Path

import mcol
import mrich
import pandas as pd
from designdb.models import (
    CataloguePriceCompoundJunctionModel,
    CataloguePriceModel,
    CompoundModel,
    CompoundTagJunctionModel,
    CompoundTagModel,
    InspirationModel,
    PoseModel,
    ReactantModel,
    ReactionModel,
    ScaffoldModel,
)
from django.db.models import Exists, OuterRef, Q
from molparse.atomtypes import formula_to_atomtype_dict
from molparse.rdkit import draw_highlighted_mol, draw_mcs
from molparse.rdkit.classify import classify_mol
from rdkit import Chem
from rdkit.Chem import Descriptors, MolFromSmarts, rdRGroupDecomposition
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcNumRings
from rdkit.Chem.Scaffolds import MurckoScaffold

from .price import Price
from .quote import Quote


class Compound:
    """A :class:`.Compound` represents a ligand/small molecule with stereochemistry
    removed and no atomic coordinates. I.e. it represents the chemical structure.
    It's name is always an InChiKey. If a compound is an elaboration it can have a
    :meth:`.Compound.scaffolds` property which is another :class:`.Compound`.
    :class:`.Compound` objects are target-agnostic and can be linked to any number of
    catalogue entries (:class:`.Quote`) or synthetic pathways (:class:`.Reaction`).

    .. attention::

            :class:`.Compound` objects should not be created directly. Instead use
            :meth:`.HIPPO.register_compound` or :meth:`.HIPPO.compounds`. See
            :doc:`getting_started` and :doc:`insert_elaborations`.

    """

    _table = 'compound'

    def __init__(self, instance: CompoundModel):
        """Compound initialisation"""

        self._instance = instance

        # caches
        self._scaffolds = None
        self._elabs = None
        self._tags = None
        self._mol = None
        self._num_heavy_atoms = None
        self._num_rings = None
        self._formula = None
        self._molecular_weight = None

    ### FACTORIES

    @classmethod
    def from_id(cls, id: int) -> 'Compound':
        """Create a :class:`.Compound` from its database ID"""
        return cls(CompoundModel.objects.get(pk=id))

    ### PROPERTIES

    @property
    def id(self) -> int:
        """Returns the compound's database ID"""
        return self._instance.pk

    @property
    def inchikey(self) -> str:
        """Returns the compound's InChiKey"""
        return self._instance.compound_inchikey

    @property
    def name(self) -> str:
        """Returns the compound's alias, or InChiKey if no alias is set"""
        if self.alias:
            return self.alias
        return self.inchikey

    @property
    def smiles(self) -> str:
        """Returns the compound's (flattened) SMILES"""
        return self._instance.compound_smiles

    @property
    def alias(self) -> str:
        """Returns the compound's alias"""
        return self._instance.compound_alias

    @alias.setter
    def alias(self, alias: str) -> None:
        """Set the compound's alias"""
        self.set_alias(alias)

    @property
    def mol(self) -> Chem.Mol | None:
        """Returns the compound's RDKit Molecule"""
        if self._mol is None:
            mol_text = self._instance.compound_mol
            if mol_text:
                self._mol = Chem.MolFromMolBlock(mol_text)
        return self._mol

    @property
    def num_heavy_atoms(self) -> int | None:
        """Get the number of heavy atoms"""
        if self._num_heavy_atoms is None and self.mol is not None:
            self._num_heavy_atoms = self.mol.GetNumHeavyAtoms()
        return self._num_heavy_atoms

    @property
    def molecular_weight(self) -> float | None:
        """Get the molecular weight"""
        if self._molecular_weight is None and self.mol is not None:
            self._molecular_weight = Descriptors.ExactMolWt(self.mol)
        return self._molecular_weight

    @property
    def num_rings(self) -> int | None:
        """Get the number of rings"""
        if self._num_rings is None and self.mol is not None:
            self._num_rings = CalcNumRings(self.mol)
        return self._num_rings

    @property
    def formula(self) -> str | None:
        """Get the chemical formula"""
        if self._formula is None and self.mol is not None:
            self._formula = CalcMolFormula(self.mol)
        return self._formula

    @property
    def atomtype_dict(self) -> dict[str, int]:
        """Get a dictionary with atomtypes as keys and corresponding
        quantities/counts as values."""

        return formula_to_atomtype_dict(self.formula)

    @property
    def num_atoms_added(self) -> int | list[int] | None:
        """Calculate the number of atoms added relative to the scaffold compound"""
        match self.num_scaffolds:
            case 0:
                mrich.error(f'{self} has no scaffold')
                return None
            case 1:
                scaffold = Compound(next(iter(self.scaffolds._queryset)))
                return self.num_heavy_atoms - scaffold.num_heavy_atoms
            case _:
                mrich.warning(f'{self} has multiple scaffolds')
                n_e = self.num_heavy_atoms
                return [
                    n_e - Compound(c).num_heavy_atoms for c in self.scaffolds._queryset
                ]

    @property
    def metadata(self) -> dict | None:
        """Returns the compound's metadata dict"""
        return self._instance.compound_metadata

    @property
    def tags(self) -> list[str]:
        """Returns the compound's tags"""
        if self._tags is None:
            self._tags = self.get_tags()
        return self._tags

    @property
    def poses(self) -> 'PoseSet':
        """Returns the compound's poses"""
        return self.get_poses()

    @property
    def best_placed_pose(self) -> 'PoseModel':
        """Returns the compound's pose with the lowest distance score"""
        return self.poses.best_placed_pose

    @property
    def num_poses(self) -> int:
        """Returns the number of associated poses"""
        return PoseModel.objects.filter(compound=self._instance).count()

    @property
    def num_reactions(self) -> int:
        """Returns the number of associated reactions (product)"""
        return ReactionModel.objects.filter(product_compound=self._instance).count()

    @property
    def num_reactant(self) -> int:
        """Returns the number of associated reactions (reactant)"""
        return ReactantModel.objects.filter(compound=self._instance).count()

    @property
    def scaffolds(self) -> 'CompoundSet | None':
        """Returns the scaffold compounds for this elaboration"""
        if self._scaffolds is None:
            ids = self.get_scaffold_ids()
            if not ids:
                return None
            from designdb.sets.compound import CompoundSet

            self._scaffolds = CompoundSet(ids, name=f'scaffolds of {self}')
        return self._scaffolds

    @property
    def num_scaffolds(self) -> int:
        """Get the number of scaffold compounds for this elaboration"""
        if scaffolds := self.scaffolds:
            return len(scaffolds)
        return 0

    @property
    def elabs(self) -> 'CompoundSet | None':
        """Returns the elaborations of this scaffold compound"""
        if self._elabs is None:
            ids = self.get_superstructure_ids()
            if not ids:
                return None
            from designdb.sets.compound import CompoundSet

            self._elabs = CompoundSet(ids, name=f'elaborations of {self}')
        return self._elabs

    @property
    def reactions(self) -> 'ReactionSet':
        """Returns the reactions resulting in this compound"""
        return self.get_reactions()

    @property
    def reaction(self) -> 'ReactionModel | None':
        """Returns the reaction resulting in this compound (warns if multiple)"""
        reactions = self.reactions
        match len(reactions):
            case 0:
                mrich.warning(f'{self} has no reactions')
                return None
            case 1:
                pass
            case _:
                mrich.warning(f'{self} has multiple reactions, returning first')
        return reactions[0]

    @property
    def dict(self) -> dict:
        """Returns a dictionary of this compound. See :meth:`.Compound.get_dict`"""
        return self.get_dict()

    @property
    def is_scaffold(self) -> bool:
        """Is this Compound the basis for any elaborations?"""
        return ScaffoldModel.objects.filter(base_compound=self._instance).exists()

    @property
    def is_elab(self) -> bool:
        """Is this Compound based on any other compound?"""
        return ScaffoldModel.objects.filter(
            superstructure_compound=self._instance
        ).exists()

    @property
    def is_product(self) -> bool:
        """Is this Compound a product of at least one reaction?"""
        return ReactionModel.objects.filter(product_compound=self._instance).exists()

    @property
    def table(self) -> str:
        """Returns the name of the database table"""
        return self._table

    ### METHODS

    def add_stock(
        self,
        amount: float,
        *,
        purity: float | None = None,
        entry: str | None = None,
        location: str | None = None,
        return_quote: bool = True,
    ) -> int | CataloguePriceModel:
        """Register a certain quantity of compound stock in the Database.

        :param amount: Amount in ``mg``
        :param purity: Purity fraction ``0 < purity <= 1``, defaults to ``None``
        :param entry: Catalogue entry identifier, defaults to ``None``
        :param location: String describing where this stock is located, defaults to
            ``None``
        :param return_quote: If ``True`` a :class:`.CataloguePriceModel` object is
            returned instead of its ID, defaults to ``True``
        :returns: The inserted :class:`.CataloguePriceModel` object or ID
        """

        assert amount

        # Find existing in-stock entries for this compound
        existing_qs = CataloguePriceModel.objects.annotate(
            has_compound=Exists(
                CataloguePriceCompoundJunctionModel.objects.filter(
                    compound=self._instance,
                    catalogue_price=OuterRef('pk'),
                )
            )
        ).filter(has_compound=True, supplier='Stock')

        # Delete entries matching this entry/purity/location
        to_delete = existing_qs.filter(
            supplier_id=entry or '',
            purity=purity,
            vendor=location or '',
        )
        deleted_count = to_delete.count()
        if deleted_count:
            CataloguePriceCompoundJunctionModel.objects.filter(
                compound=self._instance,
                catalogue_price__in=to_delete,
            ).delete()
            mrich.warning(f'Removed {deleted_count} existing In-Stock entries')

        not_deleted = existing_qs.exclude(
            supplier_id=entry or '',
            purity=purity,
            vendor=location or '',
        ).count()
        if not_deleted:
            mrich.warning(
                f'Did not remove {not_deleted} existing In-Stock entries with'
                f' differing entry/purity/location'
            )

        # Create new price entry and link to compound
        quote = CataloguePriceModel.objects.create(
            supplier='Stock',
            supplier_id=entry or '',
            vendor=location or '',
            amount=amount,
            price=0,
            currency=None,
            lead_time=0,
            purity=purity,
            catalogue_compound=None,
        )
        CataloguePriceCompoundJunctionModel.objects.create(
            compound=self._instance,
            catalogue_price=quote,
        )

        if return_quote:
            return quote
        return quote.pk

    def get_tags(self) -> list[str]:
        """Get the tags assigned to this compound"""
        return list(self._instance.tags.values_list('compound_tag_name', flat=True))

    def add_tag(self, tag: str) -> None:
        """Add a tag to this compound"""

        assert isinstance(tag, str)
        tag_obj, _ = CompoundTagModel.objects.get_or_create(compound_tag_name=tag)
        CompoundTagJunctionModel.objects.get_or_create(
            compound=self._instance, compound_tag=tag_obj
        )
        self._tags = None  # invalidate cache

    def get_quotes(
        self,
        min_amount: float | None = None,
        supplier: str | None = None,
        max_lead_time: float | None = None,
        none: str = 'quiet',
        pick_cheapest: bool = False,
        df: bool = False,
    ):
        """Get all quotes associated to this compound.
        See :meth:`.Ingredient.get_quotes`"""
        return Ingredient.get_quotes(
            compound=self._instance,
            min_amount=min_amount,
            supplier=supplier,
            max_lead_time=max_lead_time,
            none=none,
            pick_cheapest=pick_cheapest,
            df=df,
        )

    def get_reactions(
        self,
        as_reactant: bool = False,
        permitted_reactions: 'ReactionSet' = None,
        none: str = 'error',
    ) -> 'ReactionSet':
        """Get the associated :class:`.ReactionModel` objects.

        :param as_reactant: Search for reactions using this compound as a reactant,
            defaults to ``False``
        :param permitted_reactions: Filter results to this :class:`.ReactionSet`
        :param none: Unused, kept for API compatibility
        """

        from designdb.sets.reaction import ReactionSet

        if as_reactant:
            reaction_ids = list(
                ReactantModel.objects.filter(compound=self._instance).values_list(
                    'reaction_id', flat=True
                )
            )
        else:
            reaction_ids = list(
                ReactionModel.objects.filter(
                    product_compound=self._instance
                ).values_list('pk', flat=True)
            )

        if permitted_reactions:
            permitted_ids = set(permitted_reactions.ids)
            reaction_ids = [i for i in reaction_ids if i in permitted_ids]

        rset = ReactionSet(reaction_ids)
        if not as_reactant and not permitted_reactions:
            rset._name = f'reactions resulting in {str(self)}'

        return rset

    def get_poses(self) -> 'PoseSet':
        """Get the associated :class:`.PoseModel` objects."""
        from designdb.sets.pose import PoseSet

        qs = PoseModel.objects.filter(compound=self._instance)
        return PoseSet(qs, name=f"{self}'s poses")

    def get_dict(
        self,
        *,
        mol: bool = True,
        alias: bool = True,
        inchikey: bool = True,
        metadata: bool = True,
        poses: bool = True,
        num_reactant: bool = True,
        num_reactions: bool = True,
        scaffolds: bool = True,
        elabs: bool = True,
        tags: bool = True,
    ) -> dict:
        """Returns a dictionary representing this :class:`.Compound`

        :param mol: Include a ``rdkit.Chem.Mol object``, defaults to ``True``
        :param metadata: Include metadata, defaults to ``True``
        :param poses: Include IDs of associated :class:`.PoseModel` objects,
            defaults to ``True``
        :param num_reactant: include num_reactant column
        :param num_reactions: include num_reactions column
        :param scaffolds: include scaffolds column
        :param elabs: include elabs column
        :param tags: include tags column
        :returns: A dictionary
        """

        data: dict = {'id': self.id, 'smiles': self.smiles}

        if alias:
            data['alias'] = self.alias
        if inchikey:
            data['inchikey'] = self.inchikey
        if num_reactant:
            data['num_reactant'] = self.num_reactant
        if num_reactions:
            data['num_reactions'] = self.num_reactions

        if mol:
            data['mol'] = self.mol

        if scaffolds:
            data['scaffolds'] = self.scaffolds.ids if self.scaffolds else None

        if elabs:
            data['elabs'] = self.elabs.ids if self.elabs else None

        if tags:
            data['tags'] = self.tags

        if poses:
            pose_set = self.poses
            if pose_set:
                data['poses'] = pose_set.ids
                data['targets'] = pose_set.target_names

        if metadata and (metadict := self.metadata):
            for key, value in metadict.items():
                data[key] = value

        return data

    def get_recipes(
        self,
        *,
        amount: float = 1,
        debug: bool = False,
        pick_cheapest: bool = False,
        quoted_only: bool = False,
        supplier: None | str = None,
        **kwargs,
    ):
        """Get :class:`.Recipe` objects that result in this compound.
        See :meth:`.Recipe.from_compounds`"""
        from designdb.sets.compound import CompoundSet

        from .recipe import Recipe

        return Recipe.from_compounds(
            CompoundSet([self._instance.pk]),
            amount=amount,
            debug=debug,
            pick_cheapest=pick_cheapest,
            quoted_only=quoted_only,
            supplier=supplier,
            **kwargs,
        )

    def get_scaffold_ids(self) -> list[int] | None:
        """Get a list of :class:`.Compound` IDs that this object is a superstructure
        of"""
        ids = list(
            ScaffoldModel.objects.filter(
                superstructure_compound=self._instance
            ).values_list('base_compound_id', flat=True)
        )
        return ids or None

    def get_superstructure_ids(self) -> list[int] | None:
        """Get a list of :class:`.Compound` IDs that this object is a substructure of"""
        ids = list(
            ScaffoldModel.objects.filter(base_compound=self._instance).values_list(
                'superstructure_compound_id', flat=True
            )
        )
        return ids or None

    def add_scaffold(
        self, scaffold: 'Compound | CompoundModel | int', commit: bool = True
    ) -> None:
        """Add a scaffold :class:`.Compound` this molecule is derived from.

        :param scaffold: The scaffold :class:`.Compound`, :class:`.CompoundModel`,
            or its ID.
        :param commit: Unused, kept for API compatibility
        """

        if isinstance(scaffold, int):
            scaffold_model = CompoundModel.objects.get(pk=scaffold)
        elif isinstance(scaffold, CompoundModel):
            scaffold_model = scaffold
        else:
            assert scaffold._table == 'compound'
            scaffold_model = scaffold._model

        ScaffoldModel.objects.get_or_create(
            base_compound=scaffold_model,
            superstructure_compound=self._instance,
        )
        self._scaffolds = None  # invalidate cache

    def set_alias(self, alias: str, commit: bool = True) -> None:
        """Set this :class:`.Compound`'s alias.

        :param alias: The alias
        :param commit: Unused, kept for API compatibility
        """

        assert isinstance(alias, str)
        self._instance.compound_alias = alias
        self._instance.save(update_fields=['compound_alias', 'updated_on'])

    def as_ingredient(
        self,
        amount: float,
        max_lead_time: float = None,
        supplier: str = None,
        get_quote: bool = True,
        quote_none: str = 'quiet',
    ) -> 'Ingredient':
        """Convert this compound into an :class:`.Ingredient` with an associated
        amount and quote.

        :param amount: Amount in ``mg``
        :param supplier: Only search for quotes with the given supplier, defaults to
            ``None``
        :param max_lead_time: Only search for quotes with lead times less than this
            (in days), defaults to ``None``
        """

        return Ingredient.from_compound(
            compound=self._instance,
            amount=amount,
            max_lead_time=max_lead_time,
            supplier=supplier,
            get_quote=get_quote,
            quote_none=quote_none,
        )

    def draw(self, scaffolds: bool = True, align_substructure: bool = False) -> None:
        """Display this compound (and its scaffold if it has one)

        .. attention::

                This method is only intended for use within a Jupyter Notebook.

        :param align_substructure: Align the two drawings by their common
            substructure, defaults to ``False``
        """

        if scaffolds and (scaffolds := self.scaffolds):
            data = {}
            for scaffold in scaffolds:
                data[scaffold.compound_smiles] = f'C{scaffold.pk} (scaffold)'
            data[self.smiles] = str(self)

            if len(data) > 1:
                drawing = draw_mcs(
                    data,
                    align_substructure=align_substructure,
                    show_mcs=False,
                    highlight=False,
                )
                display(drawing)
            else:
                mrich.error(
                    f'Problem drawing scaffold vs {self.id=}, self referential?'
                )
                display(self.mol)
        else:
            display(self.mol)

    def draw_elabs(self) -> None:
        """Draw elaborations"""

        elabs = self.elabs

        display(self)
        display(elabs)

        if not elabs:
            mrich.error(self, 'has no elaborations')
            return self.draw()

        params = rdRGroupDecomposition.RGroupDecompositionParameters()
        params.removeAllHydrogenRGroups = False
        params.removeAllHydrogenRGroupsAndLabels = True
        params.removeHydrogensPostMatch = True

        rgd = rdRGroupDecomposition.RGroupDecomposition(
            MolFromSmarts(self.smiles), params
        )
        for c in elabs._queryset:
            mol_text = c.compound_mol
            if mol_text:
                mol = Chem.MolFromMolBlock(mol_text)
                if mol:
                    rgd.Add(mol)
        rgd.Process()

        rgroup_table = rgd.GetRGroupsAsColumns()
        core = rgroup_table['Core'][0]
        attachment_points = set()
        for rgroup in rgroup_table['Core']:
            for atom in rgroup.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    attachment_points.add(atom.GetIdx())

        drawing = draw_highlighted_mol(
            core, [(i, (0.5, 1, 0.5)) for i in attachment_points]
        )
        display(drawing)

    def classify(self, draw: bool = True) -> list[tuple[str, int]]:
        """Find RDKit Fragments within the compound molecule and draw them

        :param draw: Draw the annotated molecule, defaults to ``True``
        :returns: A list of tuples containing a descriptor (``str``) and count
            (``int``) pair
        """

        return classify_mol(self.mol, draw=draw)

    def murcko_scaffold(self, generic: bool = False) -> Chem.Mol:
        """Get the rdkit MurckoScaffold for this compound"""

        scaffold = MurckoScaffold.GetScaffoldForMol(self.mol)
        if generic:
            scaffold = MurckoScaffold.MakeScaffoldGeneric(scaffold)
        return scaffold

    def summary(
        self, metadata: bool = True, draw: bool = True, tags: bool = True
    ) -> None:
        """Print a summary of this compound

        :param metadata: Include metadata, defaults to ``True``
        :param draw: Include a 2D molecule drawing, defaults to ``True``
        """

        mrich.header(self)
        mrich.var('inchikey', self.inchikey)
        mrich.var('alias', self.alias)
        mrich.var('smiles', self.smiles)
        mrich.var('scaffolds', self.scaffolds)
        mrich.var('elabs', self.elabs)
        mrich.var('is_scaffold', self.is_scaffold)
        mrich.var('is_elab', self.is_elab)
        mrich.var('num_heavy_atoms', self.num_heavy_atoms)
        mrich.var('num_rings', self.num_rings)
        mrich.var('formula', self.formula)
        mrich.var('#reactions (product)', self.num_reactions)
        mrich.var('#reactions (reactant)', self.num_reactant)

        if tags:
            mrich.var('tags', self.tags)

        poses = self.poses
        mrich.var('#poses', len(poses))
        if poses:
            mrich.var('targets', poses.targets)

        if metadata:
            mrich.var('metadata', str(self.metadata))

        if draw:
            self.draw()

    def place(
        self,
        *,
        animal: 'HIPPO',
        reference: 'PoseModel',
        inspirations: list['PoseModel'] | None = None,
        max_ddG: float = 0.0,
        max_RMSD: float = 2.0,
        output_dir: str = 'wictor_place',
        tags: list[str] = None,
        metadata: dict = None,
        overwrite: bool = False,
    ) -> 'PoseModel | None':
        """Generate a new pose for this compound using Fragmenstein.

        :param animal: The :class:`.HIPPO` instance used to register the pose
        :param reference: Choose the :class:`.PoseModel` to use as the reference
            protein conformation
        :param inspirations: Choose the (virtual) hits to define the ligand reference,
            defaults to the ``reference``'s inspirations
        :param max_ddG: Maximum ``ddG`` value permitted, defaults to ``0.0``
        :param max_RMSD: Maximum ``RMSD`` value permitted, defaults to ``2.0``
        :param output_dir: Output directory for Fragmenstein files, defaults to
            ``wictor_place``
        :param tags: Tags to assign to the created pose, defaults to ``[]``
        :param metadata: A dictionary of metadata to assign to this compound,
            defaults to ``{}``
        :param overwrite: Delete old poses, defaults to ``False``
        """

        from fragmenstein import Wictor

        tags = tags or []
        metadata = metadata or {}

        inspirations = inspirations or reference.inspirations.all()
        target = reference.target.target_name

        inspiration_mols = []
        for insp in inspirations:
            mol_text = insp.compound.compound_mol
            if mol_text:
                mol = Chem.MolFromMolBlock(mol_text)
                if mol:
                    inspiration_mols.append(mol)

        protein_pdb_block = reference.protein_system.pdb_block_with_alt_sites

        victor = Wictor(hits=inspiration_mols, pdb_block=protein_pdb_block)
        victor.work_path = output_dir
        victor.enable_stdout(logging.CRITICAL)
        victor.place(self.smiles, long_name=self.name)

        metadata['ddG'] = (
            victor.energy_score['bound']['total_score']
            - victor.energy_score['unbound']['total_score']
        )
        metadata['RMSD'] = victor.mrmsd.mrmsd

        if metadata['ddG'] > max_ddG:
            return None
        if metadata['RMSD'] > max_RMSD:
            return None

        pose = animal.register_pose(
            compound=self,
            target=target,
            path=Path(victor.work_path) / self.name / f'{self.name}.minimised.mol',
            inspirations=inspirations,
            reference=reference,
            tags=tags,
            metadata=metadata,
        )

        if overwrite:
            PoseModel.objects.filter(compound=self._instance).exclude(
                pk=pose.pk
            ).delete()
            mrich.success(f'Successfully posed {self} (and deleted old poses)')
        else:
            mrich.success(f'Successfully posed {self}')

        return pose

    def get_inspirations(
        self, debug: bool = True, none: str = 'warning'
    ) -> 'PoseSet | None':
        """Get the fragment inspirations for this compound's poses.

        Since inspirations map :class:`.PoseModel` objects to each other, this requires
        poses to be registered for this compound.

        :returns: a :class:`.PoseSet` object
        """

        from designdb.sets.pose import PoseSet

        poses_qs = PoseModel.objects.filter(compound=self._instance)
        insp_qs = InspirationModel.objects.filter(derivative_pose__in=poses_qs)

        if not insp_qs.exists() and none in ('warning', 'warn'):
            mrich.warning('Could not determine inspirations for', self)
            return None

        derivative_ids = set(insp_qs.values_list('derivative_pose_id', flat=True))
        original_ids = set(insp_qs.values_list('original_pose_id', flat=True))

        if debug:
            mrich.debug(f'Inspirations derived from {derivative_ids}')

        inspirations = PoseSet(
            PoseModel.objects.filter(pk__in=original_ids),
            name=f'Inspirations for {self}',
        )

        return inspirations

    ### DUNDERS

    def __str__(self) -> str:
        """Unformatted string representation"""
        return f'C{self.id}'

    def __repr__(self) -> str:
        """ANSI Formatted string representation"""
        return (
            f'{mcol.bold}{mcol.underline}{self} "{self.name}"'
            f'{mcol.unbold}{mcol.ununderline}'
        )

    def __rich__(self) -> str:
        """Representation for mrich"""
        return f'[bold underline]{self} "{self.name}"'

    def __eq__(self, other) -> bool:
        """Compare compounds"""
        assert isinstance(other, Compound)
        return self.id == other.id


class Ingredient:
    """An ingredient is a :class:`.Compound` with a fixed quanitity and an attached
    quote.

    .. image:: ../images/ingredient.png
              :width: 450
              :alt: Ingredient schema

    .. attention::

            :class:`.Ingredient` objects should not be created directly. Instead use
            :meth:`.Compound.as_ingredient`.
    """

    _table = 'ingredient'

    def __init__(
        self,
        compound: CompoundModel,  # or CatalogueCompoundModel?
        amount: float,
        quote: 'Quote | CataloguePriceModel | int | None' = None,
        max_lead_time: float | None = None,
        supplier: str | None = None,
    ):
        """Ingredient initialisation

        ``quote`` may be a :class:`.Quote`, a :class:`.CataloguePriceModel`, a quote
        ID (``int``), or ``None``. When only an ID (or nothing) is stored, the quote
        is resolved lazily on first access via :attr:`.quote` -- which re-quotes the
        catalogue (estimating a price if no pack is large enough).
        """

        self._compound = compound
        self._amount = amount
        self._max_lead_time = max_lead_time
        self._supplier = supplier

        # `_quote` holds a resolved Quote component (or None); `_quote_id` holds a
        # persisted quote ID awaiting lazy resolution. Estimated quotes have no ID
        # and so are always carried as a resolved `_quote` object.
        self._quote = None
        self._quote_id = None

        if isinstance(quote, Quote):
            self._quote = quote
            self._quote_id = quote.id
        elif isinstance(quote, CataloguePriceModel):
            self._quote = Quote(quote)
            self._quote_id = quote.pk
        elif quote is not None:
            self._quote_id = int(quote)

    def __str__(self) -> str:
        """Plain string representation"""
        return f'{self.amount:.2f}mg of C{self._compound.id}'

    def __repr__(self) -> str:
        """ANSI Formatted string representation"""
        return f'{mcol.bold}{mcol.underline}{str(self)}{mcol.unbold}{mcol.ununderline}'

    def __rich__(self) -> str:
        """Representation for mrich"""
        return f'[bold underline]{str(self)}'

    def __eq__(self, other) -> bool:
        """Equality operator"""

        if self.compound != other.compound:
            return False

        return self.amount == other.amount

    def __getattr__(self, key: str):
        """For missing attributes try getting from associated :class:`.Compound`"""
        return getattr(self.compound, key)

    @classmethod
    def from_compound(
        cls,
        compound: CompoundModel,
        amount: float,
        max_lead_time: float = None,
        supplier: str = None,
        get_quote: bool = True,
        quote_none: str = 'quiet',
    ) -> 'Ingredient':
        """Convert this compound into an :class:`.Ingredient` object with an
        associated amount (in ``mg``) and :class:`.Quote` if available.

        :param amount: Amount in ``mg``
        :param supplier: Only search for quotes with the given supplier, defaults to
            ``None``
        :param max_lead_time: Only search for quotes with lead times less than this
            (in days), defaults to ``None``
        """

        if get_quote:
            # quote = self.get_quotes(
            #     pick_cheapest=True,
            #     min_amount=amount,
            #     max_lead_time=max_lead_time,
            #     supplier=supplier,
            #     none=quote_none,
            # )

            # if not quote:
            #     quote = None

            quote = cls.get_quotes(
                compound=compound,
                pick_cheapest=True,
                min_amount=amount,
                max_lead_time=max_lead_time,
                supplier=supplier,
                none=quote_none,
            )

        else:
            quote = None

        return Ingredient(
            compound=compound,
            amount=amount,
            quote=quote,
            supplier=supplier,
            max_lead_time=max_lead_time,
        )

    @classmethod
    def get_quotes(
        cls,
        compound: CompoundModel,
        min_amount: float | None = None,
        supplier: str | None = None,
        max_lead_time: float | None = None,
        none: str = 'quiet',
        pick_cheapest: bool = False,
        df: bool = False,
    ):
        """Get all quotes associated to this compound

        :param min_amount: Only return quotes with amounts greater than this,
            defaults to ``None``
        :param supplier: Only return quotes with the given supplier, defaults to
            ``None``
        :param max_lead_time: Only return quotes with lead times less than this
            (in days), defaults to ``None``
        :param none: Define the behaviour when no quotes are found. Choose `error`
            to raise print an error.
        :param pick_cheapest: If ``True`` only the cheapest :class:`.Quote` is
            returned, defaults to ``False``
        :param df: Returns a ``DataFrame`` of the quoting data, defaults to ``False``
        :returns: List of :class:`.Quote` objects, ``DataFrame``, or single
            :class:`.Quote`. See ``pick_cheapest`` and ``df`` parameters

        """

        qs = CataloguePriceModel.objects.annotate(
            has_compound=Exists(
                CataloguePriceCompoundJunctionModel.objects.filter(
                    compound=compound,
                    catalogue_price=OuterRef('pk'),
                ),
            ),
        ).filter(
            has_compound=True,
        )

        if supplier:
            if isinstance(supplier, str):
                qs = qs.filter(supplier=supplier)
            else:
                qs = qs.filter(supplier__in=supplier)

        if not qs.exists():
            return None

        if max_lead_time:
            qs = qs.filter(lead_time__lte=max_lead_time)

        estimate = None

        if min_amount:
            suitable = qs.filter(amount__gte=min_amount)

            if suitable.exists():
                qs = suitable
            else:
                # no single pack is large enough -> estimate by scaling the biggest
                # available pack's unit price (mirrors legacy Quote.combination)
                mrich.debug(
                    f'No quote available for C{compound.pk} with amount >='
                    f' {min_amount} mg. Estimating price...'
                )
                estimate = Quote.estimate(min_amount, [Quote(m) for m in qs])

        if pick_cheapest:
            if estimate is not None:
                return estimate
            cheapest = qs.order_by('price').first()
            return Quote(cheapest) if cheapest is not None else None

        if df:
            return pd.DataFrame(qs.values()).drop(columns='compound')

        return qs

    ### METHODS

    def get_cheapest_quote_id(
        self,
        min_amount: float | None = None,
        supplier: str | None = None,
        max_lead_time: float | None = None,
    ) -> int | None:
        """
        Query quotes associated to this ingredient, and return the cheapest

        :param min_amount: Only return quotes with amounts greater than this,
            defaults to ``None``
        :param supplier: Only return quotes with the given supplier, defaults to
            ``None``
        :param max_lead_time: Only return quotes with lead times less than this
            (in days), defaults to ``None``
        :param none: Define the behaviour when no quotes are found. Choose `error`
            to raise print an error.
        """

        qs = CataloguePriceModel.objects.filter(compounds=self.compound)

        if supplier:
            qs = qs.filter(supplier=supplier)

        if min_amount:
            qs = qs.filter(amount__gte=min_amount)

        if max_lead_time:
            qs = qs.filter(lead_time__lte=max_lead_time)

        return qs.order_by('price').first()

    ### PROPERTIES

    @property
    def amount(self) -> float:
        """Returns the amount (in ``mg``)"""
        return self._amount

    @property
    def id(self) -> int:
        """Returns the ID of the associated :class:`.Compound`"""
        return self._compound.id

    @property
    def compound_id(self) -> int:
        """Returns the ID of the associated :class:`.Compound`"""
        return self._compound.id

    @property
    def quote(self) -> 'Quote | None':
        """The associated :class:`.Quote`, resolved lazily.

        If a persisted quote ID is stored it is fetched; otherwise the catalogue is
        re-quoted for this ingredient's amount (estimating a price if no single pack
        is large enough). Returns ``None`` if no quote can be found or estimated.
        """
        if self._quote is None:
            if self._quote_id:
                self._quote = Quote(CataloguePriceModel.objects.get(pk=self._quote_id))
            else:
                self._quote = Ingredient.get_quotes(
                    compound=self._compound,
                    min_amount=self._amount,
                    supplier=self._supplier,
                    max_lead_time=self._max_lead_time,
                    pick_cheapest=True,
                    none='quiet',
                )
                if self._quote is not None:
                    self._quote_id = self._quote.id
        return self._quote

    @property
    def quote_id(self) -> int | None:
        """ID of the associated quote, or ``None`` (estimated/unquoted)."""
        if self._quote_id is not None:
            return self._quote_id
        quote = self.quote
        return quote.id if quote is not None else None

    @property
    def price(self) -> Price:
        """Returns the price from the associated quote, or a null Price if
        unavailable."""
        quote = self.quote
        if quote is None:
            return Price.null()
        return quote.as_price

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
        """Set the amount, invalidating the cached quote so it is re-quoted (and
        re-estimated if needed) for the new amount on next access."""

        self._amount = a
        self._quote = None
        self._quote_id = None

    @property
    def compound(self) -> CompoundModel:
        """Returns the associated :class:`.Compound`"""

        # if not self._compound:
        #     self._compound = self.db.get_compound(id=self.compound_id)
        return self._compound

    @property
    def compound_price_amount_str(self) -> str:
        """String representation including :class:`.Compound`, :class:`.Price`,
        and amount."""
        return f'{self} ({self.amount})'

    @property
    def smiles(self) -> str:
        """Returns the SMILES of the associated :class:`.Compound`"""
        return self.compound.smiles
