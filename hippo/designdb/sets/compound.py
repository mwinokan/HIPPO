import json
from collections.abc import Callable
from pathlib import Path
from statistics import mean, pvariance
from typing import TYPE_CHECKING

import mcol
import mrich
import pandas as pd
from designdb.components.compound import Compound, Ingredient
from designdb.components.reaction import Reaction
from designdb.models import (
    CompoundModel,
    CompoundTagJunctionModel,
    CompoundTagModel,
    PoseModel,
    PoseTagJunctionModel,
    ReactantModel,
    ReactionModel,
    RouteModel,
    ScaffoldModel,
)
from django.db.models import Count, Exists, OuterRef, Q
from pandas import DataFrame
from rdkit import Chem

# from rdkit.Chem import inchi
from ..utils import registration_hash_tautomer_insensitive, superparent

if TYPE_CHECKING:
    import plotly.graph_objects as go
    from designdb.sets.ingredient import IngredientSet
    from designdb.sets.pose import PoseSet
    from designdb.sets.reaction import ReactionSet
    from designdb.sets.route import RouteSet


class CompoundSet:
    """Object representing a subset of the 'compound' table in the :class:`.Database`.

    .. attention::

            :class:`.CompoundSet` objects should not be created directly. Instead use
            the :meth:`.HIPPO.compounds` property. See :doc:`getting_started` and
            :doc:`insert_elaborations`.

    Use as an iterable
    ==================

    Iterate through :class:`.CompoundModel` objects in the set:

    ::

            cset = animal.compounds[:100]

            for compound in cset:
                    ...

    Check membership
    ================

    To determine if a :class:`.CompoundModel` is present in the set:

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

    Tags and scaffold compounds can also be used to filter:

    ::

            cset = animal.compounds(tag='hits') # select compounds tagged with 'hits'
            cset = animal.compounds(scaffold=comp)  # select elaborations of comp

    """

    def __init__(
        self,
        queryset=None,
        *,
        sort: bool = True,
        name: str | None = None,
    ) -> None:
        """CompoundSet initialisation"""

        if queryset:
            if isinstance(queryset, list):
                self._queryset = CompoundModel.objects.filter(pk__in=queryset)
            else:
                self._queryset = queryset
        else:
            self._queryset = CompoundModel.objects.none()

        if sort:
            self._queryset = self._queryset.order_by('pk')

        self._name = name

    ### DUNDERS

    def __len__(self) -> int:
        """The number of compounds in this set"""
        return self._queryset.count()

    def __iter__(self):
        """Iterate through compounds in this set"""
        return iter(self._queryset)

    def __getitem__(
        self,
        key: int | slice,
    ) -> 'CompoundModel | CompoundSet':
        """Get compounds or subsets thereof from this set

        :param key: integer index or slice of indices

        """
        match key:
            case int():
                # index by position in the (ordered) set; support negative
                # indices (e.g. cset[-1] -> last compound)
                n = len(self)
                idx = key + n if key < 0 else key
                if not 0 <= idx < n:
                    raise IndexError(f'CompoundSet index out of range: {key}')
                return self._queryset[idx]

            case slice():
                # positional slice of the ordered members. Slice `.all()` (a
                # fresh, unevaluated clone) so this returns a queryset (LIMIT/
                # OFFSET) even when self._queryset is already evaluated -- an
                # evaluated queryset would otherwise slice to a list of model
                # instances. sort=False: a queryset can't be re-ordered once sliced.
                return CompoundSet(self._queryset.all()[key], sort=False)

            case _:
                raise NotImplementedError

    def __sub__(
        self,
        other: 'CompoundModel | CompoundSet | IngredientSet',
    ) -> 'CompoundSet':
        """Subtract a :class:`.CompoundModel` object or ID from this set, or subtract
        multiple at once when ``other`` is a :class:`.CompoundSet` or
        :class:`.IngredientSet`"""

        # local import to avoid the IngredientSet <-> CompoundSet cycle
        from designdb.sets.ingredient import IngredientSet

        # materialise pks so set ops stay a single flat query instead of nesting
        # pk__in subqueries, which recurses under repeated accumulation
        ids = set(self._queryset.values_list('pk', flat=True))
        match other:
            case CompoundSet() | IngredientSet():
                ids -= set(other.ids)
            case int():
                ids.discard(other)
            case _:
                raise NotImplementedError
        return CompoundSet(CompoundModel.objects.filter(pk__in=ids), sort=False)

    def __add__(
        self,
        other: 'CompoundModel | CompoundSet | IngredientSet | int',
    ) -> 'CompoundSet':
        """Add a :class:`.CompoundModel` object or ID to this set, or add multiple at
        once when ``other`` is a :class:`.CompoundSet` or :class:`.IngredientSet`"""

        # local import to avoid the IngredientSet <-> CompoundSet cycle
        from designdb.sets.ingredient import IngredientSet

        # materialise pks: keep a single flat query, avoiding nested pk__in
        # subqueries that recurse when accumulated (e.g. combining many recipes)
        ids = set(self._queryset.values_list('pk', flat=True))
        match other:
            case CompoundModel():
                ids.add(other.pk)
            case int():
                ids.add(other)
            case CompoundSet() | IngredientSet():
                ids |= set(other.ids)
            case _:
                raise NotImplementedError
        return CompoundSet(CompoundModel.objects.filter(pk__in=ids), sort=False)

    def __and__(self, other: 'CompoundSet'):
        """AND set operation, returns only compounds in both sets"""

        match other:
            case CompoundSet():
                ids = set(self._queryset.values_list('pk', flat=True)) & set(other.ids)
                return CompoundSet(CompoundModel.objects.filter(pk__in=ids), sort=False)

            case _:
                raise NotImplementedError

    def __or__(self, other: 'CompoundSet'):
        """OR set operation, returns union of both sets"""

        match other:
            case CompoundSet():
                ids = set(self._queryset.values_list('pk', flat=True)) | set(other.ids)
                return CompoundSet(CompoundModel.objects.filter(pk__in=ids), sort=False)

            case _:
                raise NotImplementedError

    def __xor__(self, other: 'CompoundSet'):
        """Exclusive OR set operation, returns all compounds in either set but
        not both"""

        match other:
            case CompoundSet():
                ids = set(self._queryset.values_list('pk', flat=True)) ^ set(other.ids)
                return CompoundSet(CompoundModel.objects.filter(pk__in=ids), sort=False)

            case _:
                raise NotImplementedError

    def __str__(self) -> str:
        """Unformatted string representation"""

        if self.name:
            s = f'{self.name}: '
        else:
            s = ''

        s += f'{{C × {len(self)}}}'

        return s

    def __repr__(self) -> str:
        """ANSI ormatted string representation"""
        return f'{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}'

    def __rich__(self) -> str:
        """Representation for mrich"""
        return f'[bold underline]{self}'

    def __contains__(self, other: 'CompoundModel | int | Ingredient'):
        """Check if a compound or ingredient is a member of this set"""
        match other:
            case CompoundModel():
                pk = other.pk
            case int():
                pk = other
            case _:
                # Ingredient (or anything exposing a compound id)
                pk = getattr(other, 'compound_id', None)
                if pk is None:
                    pk = getattr(other, 'id', None)

        if pk is None:
            return False

        return self._queryset.filter(pk=pk).exists()

    ### FILTERING

    def get_by_tag(
        self,
        tag: str,
        inverse: bool = False,
    ) -> 'CompoundSet':
        """Get all child compounds with a certain tag"""

        self._queryset = self._queryset.annotate(
            has_tag=Exists(
                CompoundTagJunctionModel.objects.filter(
                    compound=OuterRef('pk'),
                    compound_tag__compound_tag_name=tag,
                ),
            ),
        )
        if inverse:
            return CompoundSet(self._queryset.filter(has_tag=False))
        else:
            return CompoundSet(self._queryset.filter(has_tag=True))

    def get_by_metadata(self, key: str, value: str | None = None) -> 'CompoundSet':
        """Get all child compounds with by their metadata. If no value is passed, then
        simply containing the key in the metadata dictionary is sufficient

        :param key: metadata key
        :param value: metadata value (Default value = None)
        """

        q = Q(compound_metadata__has_key=key)
        if value:
            q = Q(compound_metadata__key=value)

        qs = CompoundModel.objects.filter(q)

        return CompoundSet(qs)

    def get_by_scaffold(
        self,
        scaffold: CompoundModel | int,
        none: str = 'error',
    ) -> 'CompoundSet':
        """Get all compounds that elaborate the given scaffold compound

        :param scaffold: :class:`.CompoundModel` object or ID to search by

        """

        scaffold_id = scaffold if isinstance(scaffold, int) else scaffold.id

        ids = list(
            ScaffoldModel.objects.filter(
                base_compound_id=scaffold_id,
                superstructure_compound_id__in=self.ids,
            )
            .values_list('superstructure_compound_id', flat=True)
            .distinct()
        )
        return CompoundSet(ids)

    def get_by_smiles(self, smiles: str) -> CompoundModel:
        """Get a compound in this set by SMILES, using tautomer-insensitive matching.

        :param smiles: SMILES string to search for
        :raises ValueError: if SMILES standardisation fails
        :raises CompoundModel.DoesNotExist: if no match found in this set
        """
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        try:
            sp = superparent(mol)
        except Exception as e:
            raise ValueError(f'SuperParent failed: {e}') from e

        h = registration_hash_tautomer_insensitive(sp)
        return self._queryset.get(compound_hash=h)

    def get_all_possible_reactants(
        self,
        debug: bool = False,
    ) -> 'CompoundSet':
        """Recursively searches for all the reactants that could possible be needed to
        synthesise these compounds.

        :param debug: Increased verbosity for debugging (Default value = False)

        """

        qs = CompoundModel.objects.filter(
            pk__in=ReactantModel.objects.filter(
                reaction__in=self._queryset,
            ),
        )

        seen = set(qs.values_list('id', flat=True))
        frontier = set(seen)

        while frontier:
            new = (
                set(
                    CompoundModel.objects.filter(
                        pk__in=ReactantModel.objects.filter(
                            reaction__in=self._queryset,
                        ),
                    ).values_list('pk', flat=True)
                )
                - seen
            )

            seen |= new
            frontier = new

        return CompoundSet(CompoundModel.objects.filter(pk__in=seen))

    def get_all_possible_reactions(
        self,
        debug: bool = False,
    ) -> 'ReactionSet':
        """Recursively searches for all the reactants that could possible be needed to
        synthesise these compounds.

        :param debug: Increased verbosity for debugging (Default value = False)

        """
        qs = CompoundModel.objects.filter(
            pk__in=ReactantModel.objects.filter(
                reaction__in=self._queryset,
            ),
        )

        seen = set(qs.values_list('id', flat=True))
        frontier = set(seen)

        while frontier:
            new = (
                set(
                    CompoundModel.objects.filter(
                        pk__in=ReactantModel.objects.filter(
                            reaction__in=self._queryset,
                        ),
                    ).values_list('pk', flat=True)
                )
                - seen
            )

            seen |= new
            frontier = new

        return ReactionModel.objects.filter(product__compound__in=seen)

    def get_risk_diversity(self, debug: bool = False) -> float:
        """Calculate the average spread of risk (#atoms added) for each scaffold in
        this set

        :returns: average of the standard deviations of number of atoms added for each
            scaffold

        """

        # heavy-atom count per compound, cached across scaffold edges
        nha_cache: dict[int, int | None] = {}

        def nha(compound_id: int) -> int | None:
            if compound_id not in nha_cache:
                nha_cache[compound_id] = Compound(
                    CompoundModel.objects.get(pk=compound_id)
                ).num_heavy_atoms
            return nha_cache[compound_id]

        # group the #atoms-added of each elaboration by its scaffold base
        diffs_by_base: dict[int, list[int]] = {}
        for base_id, elab_id in ScaffoldModel.objects.filter(
            superstructure_compound_id__in=self.ids
        ).values_list('base_compound_id', 'superstructure_compound_id'):
            base_nha, elab_nha = nha(base_id), nha(elab_id)
            if base_nha is None or elab_nha is None:
                continue
            diffs_by_base.setdefault(base_id, []).append(elab_nha - base_nha)

        if not diffs_by_base:
            return None

        # population variance of #atoms-added within each scaffold (0 for singletons)
        variances = [pvariance(diffs) for diffs in diffs_by_base.values()]

        if debug:
            mrich.debug(f'{variances=}')

        return mean(variances)

    def count_by_tag(
        self,
        tag: str,
    ) -> 'CompoundSet':
        """Count all child compounds with a certain tag

        :param tag: tag to filter by

        """
        return self._queryset.annotate(
            has_tag=Exists(
                CompoundTagModel.objects.filter(
                    compound=OuterRef('pk'),
                    compound_tag__compound_tag_name=tag,
                ),
            ),
        ).count()

    ### CONSOLE / NOTEBOOK OUTPUT

    def draw(self) -> None:
        """Draw a grid of all contained molecules.

        .. attention::

                This method is only intended for use within a Jupyter Notebook.

        """

        from IPython.display import display
        from molparse.rdkit import draw_grid

        compounds = [Compound(c) for c in self._queryset]
        mols = [c.mol for c in compounds]
        labels = [str(c) for c in compounds]

        display(draw_grid(mols, labels=labels))

    def grid(self) -> None:
        """Draw a grid of all contained molecules.

        .. attention::

                This method is only intended for use within a Jupyter Notebook.

        """

        self.draw()

    def summary(self, return_df: bool = False) -> None:
        """Print a summary of this compound set"""

        mrich.header(self)

        ids = list(self.ids)

        # compound-tag counts: distinct compounds per compound tag
        comp_rows = [
            dict(tag=name, num_compounds=n)
            for name, n in CompoundTagJunctionModel.objects.filter(compound_id__in=ids)
            .values_list('compound_tag__compound_tag_name')
            .annotate(n=Count('compound', distinct=True))
            .order_by('compound_tag__compound_tag_name')
        ]
        df = DataFrame(comp_rows)
        if len(df):
            df = df.set_index('tag')
        else:
            df = DataFrame(columns=['num_compounds'])
            df.index.name = 'tag'

        # pose-tag counts over the poses of these compounds
        for name, num_poses, num_posed in (
            PoseTagJunctionModel.objects.filter(pose__compound_id__in=ids)
            .values_list('pose_tag__pose_tag_name')
            .annotate(
                num_poses=Count('pose', distinct=True),
                num_posed=Count('pose__compound', distinct=True),
            )
            .order_by('pose_tag__pose_tag_name')
        ):
            df.loc[name, 'num_poses'] = num_poses
            df.loc[name, 'num_posed_compounds'] = num_posed

        df.loc['TOTAL', 'num_compounds'] = len(self)
        df.loc['TOTAL', 'num_poses'] = self.num_poses
        df.loc['TOTAL', 'num_posed_compounds'] = (
            PoseModel.objects.filter(compound_id__in=ids)
            .values('compound_id')
            .distinct()
            .count()
        )

        df = df.fillna(0)
        df = df.astype(int)

        if return_df:
            return df
        else:
            mrich.print(df)

    def interactive(
        self,
        function: Callable | None = None,
    ) -> None:
        """Creates a ipywidget to interactively navigate this PoseSet."""

        from IPython.display import display
        from ipywidgets import (
            BoundedIntText,
            Checkbox,
            GridBox,
            Layout,
            VBox,
            interactive,
            interactive_output,
        )

        if function:

            def widget(i):
                """interactive function widget"""
                compound = self[i]
                display(compound)
                function(compound)

            return interactive(
                widget,
                i=BoundedIntText(
                    value=0,
                    min=0,
                    max=len(self) - 1,
                    step=1,
                    description=f'Comp (/{len(self)}):',
                    disabled=False,
                ),
            )

        else:
            a = BoundedIntText(
                value=0,
                min=0,
                max=len(self) - 1,
                step=1,
                description=f'Comp (/{len(self)}):',
                disabled=False,
            )

            b = Checkbox(description='Name', value=True)
            c = Checkbox(description='Summary', value=False)
            d = Checkbox(description='2D', value=True)
            e = Checkbox(description='Poses', value=False)
            f = Checkbox(description='Reactions', value=False)
            g = Checkbox(description='Tags', value=False)
            h = Checkbox(description='Quotes', value=False)
            i = Checkbox(description='Metadata', value=False)
            j = Checkbox(description='Classify', value=False)

            ui1 = GridBox(
                [b, c, d], layout=Layout(grid_template_columns='repeat(3, 100px)')
            )
            ui2 = GridBox(
                [e, f, g], layout=Layout(grid_template_columns='repeat(3, 100px)')
            )
            ui3 = GridBox(
                [h, i, j], layout=Layout(grid_template_columns='repeat(3, 100px)')
            )
            ui = VBox([a, ui1, ui2, ui3])

            def widget(
                i,
                name: bool = True,
                summary: bool = True,
                draw: bool = True,
                poses: bool = True,
                reactions: bool = True,
                tags: bool = True,
                quotes: bool = True,
                metadata: bool = True,
                classify: bool = True,
            ):
                """interactive default widget"""
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
                    comp.summary(metadata=False, draw=False, tags=False)

                if draw:
                    comp.draw()

                if poses and (pset := comp.poses):
                    for p in pset:
                        mrich.print(p)
                    pset.draw()

                if reactions and (reactions := comp.reactions):
                    for r in reactions:
                        mrich.print(r)
                        r.draw()

                if tags:
                    mrich.title('Tags')
                    mrich.print(comp.tags)

                if quotes:
                    mrich.title('Quotes')
                    display(comp.get_quotes(df=True))

                if metadata:
                    mrich.title('Metadata:')
                    mrich.print(comp.metadata)

                if classify:
                    mrich.title('Classification:')
                    comp.classify()

            out = interactive_output(
                widget,
                {
                    'i': a,
                    'name': b,
                    'summary': c,
                    'draw': d,
                    'poses': e,
                    'reactions': f,
                    'tags': g,
                    'quotes': h,
                    'metadata': i,
                    'classify': j,
                },
            )

            display(ui, out)

    def tag_summary(self) -> 'pd.DataFrame':
        """Print a summary table of tags with compound counts"""

        data = [
            dict(tag=name, num_compounds=n)
            for name, n in CompoundTagJunctionModel.objects.filter(
                compound_id__in=list(self.ids)
            )
            .values_list('compound_tag__compound_tag_name')
            .annotate(n=Count('compound', distinct=True))
            .order_by('compound_tag__compound_tag_name')
        ]

        df = DataFrame(data)
        if len(df):
            df = df.set_index('tag').astype(int)

        mrich.print(df)

        return df

    ### OTHER METHODS

    def get_recipes(
        self,
        amount: float = 1,
        debug: bool = False,
        pick_cheapest: bool = False,
        permitted_reactions: 'ReactionSet | None' = None,
        quoted_only: bool = False,
        supplier: None | str = None,
        **kwargs,
    ):
        """Generate the :class:`.Recipe` to make these compounds.

        See :meth:`.RecipeService.from_compounds`
        """

        # convenience bridge to the service layer
        from designdb.services.recipe import RecipeService

        return RecipeService.from_compounds(
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
        permitted_reactions: 'None | ReactionSet' = None,
        return_ids: bool = False,
        debug: bool = True,
    ) -> 'RouteSet':
        """Get a RoutSet to products in this set.

        :param permitted_reactions: optionally restrict reactions to those in this
            :class:`.ReactionSet`

        """

        from designdb.models import ComponentModel, RouteModel

        from .route import RouteSet

        base_qs = RouteModel.objects.filter(product_compound_id__in=list(self.ids))

        if permitted_reactions is not None:
            permitted = set(permitted_reactions.ids)

            if debug:
                mrich.debug('Querying database for routes')

            # reaction components (component_type == 1) grouped per route
            rows = ComponentModel.objects.filter(
                route__in=base_qs,
                component_type=1,
            ).values_list('route_id', 'component_ref')

            route_reactions: dict[int, set[int]] = {}
            for route_id, reaction_id in rows:
                route_reactions.setdefault(route_id, set()).add(reaction_id)

            if debug:
                mrich.debug('Checking availability')

            available_routes = [
                route_id
                for route_id, reactions in route_reactions.items()
                if reactions <= permitted
            ]

            if return_ids:
                return available_routes

            return RouteSet.from_ids(available_routes)

        route_ids = list(base_qs.values_list('id', flat=True))

        if return_ids:
            return route_ids

        return RouteSet.from_ids(route_ids)

    def copy(self) -> 'CompoundSet':
        """Returns a copy of this set"""
        return CompoundSet(self.ids)

    def get_df(
        self,
        smiles: bool = True,
        inchikey: bool = False,
        alias: bool = True,
        mol: bool = False,
        metadata: bool = False,
        expand_metadata: bool = True,
        poses: bool = False,
        num_reactant: bool = False,
        num_reactions: bool = False,
        num_poses: bool = False,
        tags: bool = False,
        scaffolds: bool = False,
        elabs: bool = False,
        routes: bool = False,
        debug: bool = False,
        **kwargs,
    ) -> 'DataFrame':
        """Get a DataFrame representation of this set

        :param smiles: include SMILES column (Default value = True)
        :param inchikey: include InChIKey column (Default value = False)
        :param alias: include alias column (Default value = True)
        :param mol: include ``rdkit.Chem.Mol`` in output (Default value = False)
        :param metadata: include metadata in output (Default value = False)
        :param expand_metadata: create separate column for each metadata key
            (Default value = True)
        :param poses: include poses in output (Default value = False)
        :param num_reactant: include num_poses column
        :param num_reactant: include num_reactant column (number of reactions where
            compound is a reactant)
        :param num_reactions: include num_reactions column (number of reactions where
            compound is a product)
        :param tags: include tags column
        :param scaffolds: include scaffolds column
        :param elabs: include elabs column

        """

        if mol:
            raise NotImplementedError(
                'get_df(mol=True) requires the RDKit cartridge (mol_to_binary_mol)'
            )

        ids = list(self.ids)

        # base columns straight off the CompoundModel rows
        fields = ['id']
        if smiles:
            fields.append('compound_smiles')
        if inchikey:
            fields.append('compound_inchikey')
        if alias:
            fields.append('compound_alias')
        if metadata:
            fields.append('compound_metadata')

        if debug:
            mrich.debug('querying...')

        rows = {
            r['id']: r for r in CompoundModel.objects.filter(pk__in=ids).values(*fields)
        }

        data = []
        for cid in ids:
            row = rows.get(cid)
            if row is None:
                continue

            d = dict(id=cid)

            if smiles:
                d['smiles'] = row['compound_smiles']

            if inchikey:
                d['inchikey'] = row['compound_inchikey']

            if alias:
                d['alias'] = row['compound_alias']

            # compound_metadata is JSON stored in a TextField
            if metadata and (meta_str := row['compound_metadata']):
                meta_dict = json.loads(meta_str)

                if expand_metadata:
                    for k, v in meta_dict.items():
                        d[k] = v
                else:
                    d['metadata'] = meta_dict

            data.append(d)

        df = DataFrame(data)

        if not data:
            return df

        if poses or num_poses:
            if debug:
                mrich.debug('adding pose column')

            lookup: dict[int, set] = {}
            for cid, pid in PoseModel.objects.filter(compound_id__in=ids).values_list(
                'compound_id', 'id'
            ):
                lookup.setdefault(cid, set()).add(pid)

            if poses:
                df['poses'] = df['id'].apply(lambda x: lookup.get(x, set()))
            if num_poses:
                df['num_poses'] = df['id'].apply(lambda x: len(lookup.get(x, set())))

        if num_reactant:
            if debug:
                mrich.debug('adding num_reactant column')
            counts: dict[int, int] = {}
            for cid in ReactantModel.objects.filter(compound_id__in=ids).values_list(
                'compound_id', flat=True
            ):
                counts[cid] = counts.get(cid, 0) + 1
            df['num_reactant'] = df['id'].apply(lambda x: counts.get(x, 0))

        if num_reactions:
            if debug:
                mrich.debug('adding num_reactions column')
            counts = {}
            for cid in ReactionModel.objects.filter(
                product_compound_id__in=ids
            ).values_list('product_compound_id', flat=True):
                counts[cid] = counts.get(cid, 0) + 1
            df['num_reactions'] = df['id'].apply(lambda x: counts.get(x, 0))

        if scaffolds:
            if debug:
                mrich.debug('adding scaffolds column')
            lookup = {}
            for sup_id, base_id in ScaffoldModel.objects.filter(
                superstructure_compound_id__in=ids
            ).values_list('superstructure_compound_id', 'base_compound_id'):
                lookup.setdefault(sup_id, set()).add(base_id)
            df['scaffolds'] = df['id'].apply(lambda x: lookup.get(x, set()))

        if elabs:
            if debug:
                mrich.debug('adding elabs column')
            lookup = {}
            for base_id, sup_id in ScaffoldModel.objects.filter(
                base_compound_id__in=ids
            ).values_list('base_compound_id', 'superstructure_compound_id'):
                lookup.setdefault(base_id, set()).add(sup_id)
            df['elabs'] = df['id'].apply(lambda x: lookup.get(x, set()))

        if tags:
            if debug:
                mrich.debug('adding tag column')
            lookup = {}
            for cid, name in CompoundTagJunctionModel.objects.filter(
                compound_id__in=ids
            ).values_list('compound_id', 'compound_tag__compound_tag_name'):
                lookup.setdefault(cid, set()).add(name)
            df['tags'] = df['id'].apply(lambda x: lookup.get(x, set()))

        if routes:
            if debug:
                mrich.debug('adding route column')
            lookup = {}
            for pid, rid in RouteModel.objects.filter(
                product_compound_id__in=ids
            ).values_list('product_compound_id', 'id'):
                lookup.setdefault(pid, set()).add(rid)
            df['routes'] = df['id'].apply(lambda x: lookup.get(x, set()))

        df = df.set_index('id')

        return df

    def get_quoted(
        self,
        *,
        supplier: str = 'any',
    ) -> 'CompoundSet':
        """Get all member compounds that have a catalogue quote.

        Quotes live in ``catalogue_prices`` and are linked to compounds via the
        ``compound_catalogue_map`` junction (see :class:`.QuoteService`).

        :param supplier: restrict to this supplier, or ``'any'`` (default)
        """

        from designdb.models import CataloguePriceCompoundJunctionModel

        qs = CataloguePriceCompoundJunctionModel.objects.filter(
            compound_id__in=list(self.ids)
        )
        if supplier != 'any':
            qs = qs.filter(catalogue_price__supplier=supplier)

        quoted = set(qs.values_list('compound_id', flat=True).distinct())
        return CompoundSet([i for i in self.ids if i in quoted])

    def get_unquoted(
        self,
        *,
        supplier: str = 'any',
    ) -> 'CompoundSet':
        """Get all member compounds that do not have a quote from given supplier

        :param supplier: supplier name (Default value = 'any')

        """

        quoted = self.get_quoted(supplier=supplier)
        return self - quoted

    def get_dict(self) -> dict:
        """Get a dictionary object with all serialisable data needed to reconstruct
        this set"""
        return dict(indices=list(self.indices))

    def write_smiles_csv(
        self, file: str, tags: bool = True, split_tags: bool = True
    ) -> None:
        """Write a CSV of the smiles contained in this set to a file

        :param file: path of the CSV file
        :param tags: include tags in output
        :param split_tags: split tags into separate columns

        """
        df = self.get_df(smiles=True, alias=False, tags=tags).reset_index()

        # explode the per-compound tag sets into a boolean column per tag
        if tags and split_tags and 'tags' in df.columns:
            all_tags: set[str] = set().union(*df['tags']) if len(df) else set()
            for tag in sorted(all_tags):
                df[tag] = df['tags'].apply(lambda s, t=tag: t in s)
            df = df.drop(columns=['tags'])

        mrich.writing(file)
        df.to_csv(file, index=False)

    def write_postera_csv(
        self,
        file,
        *,
        supplier: str = 'Enamine',
        prefix: str = 'fragment',
    ) -> None:
        """Write a CSV formatted for upload to Postera's Manifold

        :param file: path of the CSV file
        :param supplier: supplier to use for quotes, (Default value = 'Enamine')
        :param prefix: prefix to metadata columns, (Default value = 'fragment')

        """

        from datetime import date as dt

        from pandas import DataFrame

        if prefix:
            prefix = f'{prefix}_'

        data = []

        for c in mrich.track(self, prefix='Creating DataFrame'):
            # get props
            smiles = c.smiles
            tags = c.tags
            metadata = c.metadata
            poses = c.poses
            scaffold = c.scaffold

            # method
            assert len(tags) == 1, c
            method = tags[0]

            # date
            date = dt.today()

            # author
            assert 'author' in metadata, c
            author = metadata['author']

            match len(poses):
                case 1:
                    pose = poses[0]
                case 0:
                    mrich.warning(f'{c} has no poses')
                    assert scaffold
                    pose = scaffold.poses[0]
                case _:
                    mrich.warning(f'{c} has multiple poses')
                    pose = poses[0]

            # extract inspirations
            inspirations = pose.inspirations
            inspiration_names = ','.join(inspirations.names)
            inspiration_smiles = '.'.join(inspirations.smiles)

            # quote info
            quotes = c.get_quotes(supplier=supplier)
            assert len(quotes) == 1, c
            quote = quotes[0]
            catalog_id = quote.entry
            catalog_price = quote.price
            catalog_lead_time = quote.lead_time

            # hippo string
            hippo_str = f'compound={c.id}, pose={pose.id}'

            # create row
            data.append(
                {
                    'SMILES': smiles,
                    f'{prefix}HIPPO_IDs': hippo_str,
                    f'{prefix}method': method,
                    f'{prefix}export_date': date,
                    f'{prefix}author': author,
                    f'{prefix}inspiration_names': inspiration_names,
                    f'{prefix}inspiration_SMILES': inspiration_smiles,
                    f'{prefix}supplier': supplier,
                    f'{prefix}supplier_catalogue': quote.catalogue,
                    f'{prefix}supplier_ID': catalog_id,
                    f'{prefix}supplier_price': catalog_price,
                    f'{prefix}supplier_lead_time': catalog_lead_time,
                }
            )

        df = DataFrame(data)

        mrich.writing(file)
        df.to_csv(file, index=False)

        return df

    def write_CAR_csv(
        self,
        file: 'str | Path',
        amount: float = 1,  # in mg
        return_df: bool = False,
        # pick_cheapest: bool = False,
        quoted_only: bool = False,
        get_ingredient_quotes: bool = True,
        **kwargs,
    ) -> 'DataFrame | None':
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

        # avoiding circular imports
        from designdb.recipe import Recipe

        file = str(Path(file).resolve())

        rows = []

        for r_id in mrich.track(self.reaction_ids, prefix='Solving compound recipes'):
            reaction = ReactionModel.objects.get(pk=r_id)

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
                    'target-names': str(product.compound),
                    'no-steps': 0,
                    'concentration-required-mM': None,
                    'amount-required-uL': None,
                    'batch-tag': None,
                }

                for i, reaction_model in enumerate(sub_recipe.reactions):
                    reaction = Reaction(reaction_model)
                    i = i + 1

                    row['no-steps'] += 1

                    match len(reaction.reactants):
                        case 1:
                            row[f'reactant-1-{i}'] = reaction.reactants[0].smiles
                            row[f'reactant-2-{i}'] = None
                        case 2:
                            row[f'reactant-1-{i}'] = reaction.reactants[0].smiles
                            row[f'reactant-2-{i}'] = reaction.reactants[1].smiles
                        case _:
                            raise NotImplementedError(
                                f'Unsupported number of reactants for'
                                f' {reaction=}: {len(reaction.reactants)}'
                            )

                    row[f'reaction-product-smiles-{i}'] = reaction.product.smiles
                    row[f'reaction-name-{i}'] = reaction.type
                    row[f'reaction-recipe-{i}'] = None
                    row[f'reaction-groupby-column-{i}'] = None
                    # row[f'reaction-id-{i}'] = int(reaction.id)

                rows.append(row)

        df = DataFrame(rows)

        df = df.convert_dtypes()

        for n_steps in set(df['no-steps']):
            subset = df[df['no-steps'] == n_steps]
            this_file = file.replace('.csv', f'_{n_steps}steps.csv')
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

        compound_tag, _ = CompoundTagModel.objects.get_or_create(compound_tag_name=tag)

        CompoundTagJunctionModel.objects.bulk_create(
            [
                CompoundTagJunctionModel(compound=compound, compound_tag=compound_tag)
                for compound in self._queryset
            ],
            ignore_conflicts=True,
        )

        mrich.print(f'Tagged {self} w/ "{tag}"')

    def plot_tsnee(self, **kwargs) -> 'go.Figure':
        """Plot a tanimoto similarity plot of these compounds"""
        from .plotting import plot_compound_tsnee

        return plot_compound_tsnee(self, **kwargs)

    def as_ingredientset(
        self,
        amount: float | list[float] = 1,
        supplier: str | list | None = None,
    ) -> 'IngredientSet':
        """Get an :class:`.IngredientSet` for these compounds"""
        from designdb.sets.ingredient import IngredientSet

        return IngredientSet.from_compounds(
            compounds=self, amount=amount, supplier=supplier
        )

    def split_by_scaffolds(self) -> 'dict[CompoundSet, CompoundSet]':
        """Split this set into subsets clustered by scaffold compound

        Maps a single-compound :class:`.CompoundSet` for each scaffold base to a
        :class:`.CompoundSet` of its elaborations within this set.
        """
        clusters: dict[int, set] = {}
        for base_id, sup_id in ScaffoldModel.objects.filter(
            superstructure_compound_id__in=self.ids
        ).values_list('base_compound_id', 'superstructure_compound_id'):
            clusters.setdefault(base_id, set()).add(sup_id)

        return {
            CompoundSet([base_id]): CompoundSet(list(elab_ids))
            for base_id, elab_ids in clusters.items()
        }

    def despaghettify(
        self,
        register_missing_routes: bool = True,
        supplier='Enamine',
    ) -> 'CompoundSet':
        """Reduce this set to only compounds that elaborate a single reactant at a time.
        Requires routes to be present in the database."""

        if register_missing_routes:
            mrich.debug('registering_missing_routes...')
            route_lookup = self.register_missing_routes(
                missing_only=True, supplier=supplier
            )

        mrich.debug('clustering by scaffold...')
        clustered = self.split_by_scaffolds()

        n = len(clustered)
        mrich.var('#clusters', n)

        mrich.debug('getting route lookup...')
        route_lookup = self.db.get_product_id_routes_dict()

        mrich.debug('getting reactant lookup...')
        reactant_lookup = self.db.get_route_id_reactant_ids_dict()

        keep = set()
        for i, (cluster, elabs) in enumerate(clustered.items()):
            for scaffold in cluster:
                mrich.debug(
                    f'{i}/{n}',
                    'scaffold:',
                    scaffold.id,
                    '#elabs:',
                    len(elabs),
                    '#kept:',
                    len(keep),
                )

                route_ids = route_lookup.get(scaffold.id)

                if not route_ids:
                    mrich.error(f'scaffold {scaffold} has no routes')
                    continue

                elif len(route_ids) > 1:
                    mrich.warning(f'scaffold {scaffold} has multiple routes')

                for route_id in route_ids:
                    scaffold_reactants = reactant_lookup[route_id]

                    for elab in elabs:
                        route_ids = route_lookup.get(elab.id, set())

                        if len(route_ids) != 1:
                            mrich.error(f'elab {elab.id} has {route_ids=}')
                            continue

                        reactants = reactant_lookup[list(route_ids)[0]]

                        common = scaffold_reactants & reactants

                        if len(common) == len(scaffold_reactants) - 1:
                            keep.add(elab.id)

        return CompoundSet(self.db, keep)

    def register_missing_routes(
        self, missing_only: bool = True, supplier: str = 'Enamine'
    ) -> None:
        """Calculate missing routes to compounds in this set"""

        if missing_only:
            records = self.db.select_where(
                table='route',
                key=f'route_product IN {self.str_ids}',
                query='route_product',
                multiple=True,
            )
            existing = set(i for (i,) in records)
            missing = set(self.ids) - existing
            return CompoundSet(self.db, missing).register_missing_routes(
                missing_only=False, supplier=supplier
            )

        mrich.var('#compounds', len(self))

        for c in mrich.track(self, total=len(self)):
            try:
                reactions = c.reactions
            except Exception as e:
                mrich.error(f"Error getting {c}'s reactions", e)
                continue

            for reaction in reactions:
                try:
                    recipes = reaction.get_recipes(supplier=supplier)
                except Exception as e:
                    mrich.error(f"Error getting {reaction}'s ({c}) recipes", e)
                    continue

                for recipe in recipes:
                    route = self.db.register_route(recipe=recipe)

                    mrich.print(f'registered {route=}')

        self.db.prune_duplicate_routes()

    ### PROPERTIES

    @property
    def queryset(self):
        """Associated :class:`.Database` object"""
        return self._queryset

    @property
    def indices(self) -> list[int]:
        """Returns the ids of compounds in this set"""
        return self._queryset.values_list('id', flat=True)

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
        return list(self._queryset.values_list('compound_alias', flat=True))

    @property
    def smiles(self) -> list[str]:
        """Returns the smiles of child compounds"""
        return list(self._queryset.values_list('compound_smiles', flat=True))

    @property
    def mols(self) -> 'list[Chem.Mol]':
        """Returns the RDKit molecules of child compounds

        Uses each compound's stored MolBlock (``compound_mol``) where available,
        falling back to parsing ``compound_smiles``.
        """
        mols = []
        for molblock, smiles in self._queryset.values_list(
            'compound_mol', 'compound_smiles'
        ):
            if molblock:
                mols.append(Chem.MolFromMolBlock(molblock))
            elif smiles:
                mols.append(Chem.MolFromSmiles(smiles))
            else:
                mols.append(None)
        return mols

    @property
    def inchikeys(self) -> list[str]:
        """Returns the inchikeys of compounds in this set"""
        return list(self._queryset.values_list('compound_inchikey', flat=True))

    @property
    def tags(self) -> set[str]:
        """Returns the set of unique tags present in this compound set"""
        return set(
            CompoundTagJunctionModel.objects.filter(compound_id__in=self.ids)
            .values_list('compound_tag__compound_tag_name', flat=True)
            .distinct()
        )

    @property
    def num_poses(self) -> int:
        """Count the poses associated to this set of compounds"""
        return PoseModel.objects.filter(compound_id__in=self.ids).count()

    @property
    def poses(self) -> 'PoseSet':
        """Get the poses associated to this set of compounds"""
        from .pose import PoseSet

        return PoseSet(PoseModel.objects.filter(compound_id__in=self.ids))

    def best_placed_poses(
        self,
        scoring_method: str,
        version: str | None = None,
        inverse: bool = False,
    ) -> 'PoseSet':
        """Get the best-scoring pose for each compound in this set.

        Delegates to :meth:`.PoseSet.get_best_scoring_poses_per_compound`.

        :param scoring_method: ``ScoringMethodModel.method_name`` to rank by
        :param version: ``ScoringMethodModel.method_version`` — required when
            multiple versions of the same method exist
        :param inverse: if ``True``, higher score is better (default: lower is better)
        """
        return self.poses.get_best_scoring_poses_per_compound(
            scoring_method, version=version, inverse=inverse
        )

    @property
    def str_ids(self) -> str:
        """Return an SQL formatted tuple string of the :class:`.CompoundModel` IDs"""
        return str(tuple(self.ids)).replace(',)', ')')

    @property
    def num_heavy_atoms(self) -> int:
        """Get the total number of heavy atoms"""
        return sum(
            n for c in self._queryset if (n := Compound(c).num_heavy_atoms) is not None
        )

    @property
    def num_rings(self):
        """Get the total number of molecular rings"""
        return sum(
            n for c in self._queryset if (n := Compound(c).num_rings) is not None
        )

    @property
    def formula(self) -> str:
        """Get the combined chemical formula for all compounds"""
        from molparse.atomtypes import atomtype_dict_to_formula

        return atomtype_dict_to_formula(self.atomtype_dict)

    @property
    def atomtype_dict(self) -> dict[str, int]:
        """Get a dictionary with atomtypes as keys and corresponding
        quantities/counts as values"""
        from molparse.atomtypes import combine_atomtype_dicts

        atomtype_dicts = [Compound(c).atomtype_dict for c in self._queryset]
        return combine_atomtype_dicts(atomtype_dicts)

    @property
    def num_atoms_added(self) -> list:
        """Calculate the number of heavy atoms added w.r.t the scaffold(s)

        Delegates to :attr:`.Compound.num_atoms_added` per member, so each entry
        is an ``int`` (single scaffold), a list of ints (multiple scaffolds), or
        ``None`` (no scaffold), aligned with this set's compounds.
        """
        return [Compound(c).num_atoms_added for c in self._queryset]

    @property
    def avg_num_atoms_added(self) -> float:
        """Calculate the average number of heavy atoms added w.r.t the scaffold

        :returns: mean number of atoms added across members which have a scaffold
            (``0.0`` if none do)
        """
        values: list[int] = []
        for v in self.num_atoms_added:
            if v is None:
                continue
            if isinstance(v, list):
                values.extend(v)
            else:
                values.append(v)
        return mean(values) if values else 0.0

    @property
    def risk_diversity(self) -> float:
        """Calculate the average spread of risk (#atoms added) for each scaffold in
        this set

        :returns: average of the standard deviations of number of atoms added for each
            scaffold

        """

        return self.get_risk_diversity()

    @property
    def elaboration_balance(self) -> float:
        """Measure of how evenly elaborations are distributed across scaffolds in
        this set"""

        # number of elaborations (in this set) per scaffold base compound
        counts_by_base: dict[int, int] = {}
        for base_id in ScaffoldModel.objects.filter(
            superstructure_compound_id__in=self.ids
        ).values_list('base_compound_id', flat=True):
            counts_by_base[base_id] = counts_by_base.get(base_id, 0) + 1

        # optional dependency, isolated so module import never depends on it
        from hirsch import hirsch

        return hirsch(list(counts_by_base.values()))

    @property
    def num_scaffolds_elaborated(self) -> int:
        """Count the number of scaffold compounds that have at least one elaboration in
        this set

        :returns: number of scaffold compounds

        """
        return len(self.scaffold_ids)

    @property
    def scaffolds(self) -> 'CompoundSet':
        """Get the scaffold compounds that have at least one elaboration in this set

        :returns: :class:`.CompoundSet`

        """
        return CompoundSet(self.scaffold_ids)

    @property
    def scaffold_ids(self) -> list[int]:
        """Return a list of :class:`.CompoundModel` ID's for scaffolds of this set"""
        return list(
            ScaffoldModel.objects.filter(superstructure_compound_id__in=self.ids)
            .values_list('base_compound_id', flat=True)
            .distinct()
        )

    @property
    def num_scaffolds(self) -> int:
        """Return a count of scaffolds of this set"""
        return len(self.scaffold_ids)

    @property
    def elabs(self) -> 'CompoundSet':
        """Returns a :class:`.CompoundSet` of all compounds that are an elaboration
        of a scaffold in this set"""

        ids = list(
            ScaffoldModel.objects.filter(
                base_compound_id__in=self.ids,
                superstructure_compound_id__isnull=False,
            )
            .values_list('superstructure_compound_id', flat=True)
            .distinct()
        )
        return CompoundSet(ids)

    @property
    def num_elabs(self) -> int:
        """Return a count of elaborations of this set"""
        return len(self.elabs)

    @property
    def elab_df(self) -> 'pd.DataFrame':
        """Get a DataFrame summarising the elaborations in this CompoundSet.

        Groups this set's compounds by the scaffold(s) they elaborate; a compound
        with multiple scaffolds is listed under each.
        """
        # base scaffold compound -> elaboration compound ids within this set
        clusters: dict[int, set] = {}
        for base_id, sup_id in ScaffoldModel.objects.filter(
            superstructure_compound_id__in=self.ids
        ).values_list('base_compound_id', 'superstructure_compound_id'):
            clusters.setdefault(base_id, set()).add(sup_id)

        data = []
        for base_id, elab_ids in clusters.items():
            elabs = CompoundSet(list(elab_ids))
            data.append(
                dict(
                    scaffold_id=base_id,
                    scaffold_compound=Compound(CompoundModel.objects.get(pk=base_id)),
                    elabs=elabs,
                    num_elabs=len(elabs),
                )
            )

        return DataFrame(data)

    @property
    def id_num_poses_dict(self) -> dict[int, int]:
        """Get a dictionary mapping compound ids to the number of poses"""

        counts = dict(
            PoseModel.objects.filter(compound_id__in=self.ids)
            .values_list('compound_id')
            .annotate(n=Count('id'))
        )

        return {cid: counts.get(cid, 0) for cid in self.ids}

    @property
    def reaction_ids(self) -> list[int]:
        """Returns a list of :class:`.ReactionModel` IDs that result in members of
        this set"""
        return list(
            ReactionModel.objects.filter(product_compound_id__in=self.ids)
            .values_list('id', flat=True)
            .distinct()
        )
