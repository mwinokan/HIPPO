import json
from collections.abc import Callable
from pathlib import Path

import mcol
import mrich
import pandas as pd
from designdb.components.compound import Ingredient
from designdb.components.price import Price
from designdb.models import (
    CataloguePriceModel,
    CompoundModel,
    CompoundTagJunctionModel,
    CompoundTagModel,
    ReactantModel,
    ReactionModel,
)
from django.db.models import Exists, OuterRef, Q
from pandas import DataFrame, concat, isna
from rdkit import Chem
# from rdkit.Chem import inchi
from rdkit.Chem import Mol

from ..utils import registration_hash_tautomer_insensitive, superparent


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
                index = self.indices[key]
                try:
                    return CompoundModel.objects.get(id=index)
                except CompoundModel.DoesNotExist:
                    raise CompoundModel.DoesNotExist from exc

            case slice():
                return CompoundSet(CompoundModel.objects.filter(pk__in=key))

            case _:
                raise NotImplementedError

    def __sub__(
        self,
        other: 'CompoundModel | CompoundSet | IngredientSet',
    ) -> 'CompoundSet':
        """Subtract a :class:`.CompoundModel` object or ID from this set, or subtract
        multiple at once when ``other`` is a :class:`.CompoundSet` or
        :class:`.IngredientSet`"""

        match other:
            case CompoundSet():
                return CompoundSet(
                    CompoundModel.objects.filter(
                        Q(pk__in=self._queryset) & ~Q(pk__in=other.queryset)
                    ),
                    sort=False,
                )
            case int():
                return CompoundSet(
                    CompoundModel.objects.filter(
                        Q(pk__in=self._queryset) & ~Q(pk=other.pk)
                    ),
                    sort=False,
                )

    def __add__(
        self,
        other: 'CompoundModel | CompoundSet | IngredientSet | int',
    ) -> 'CompoundSet':
        """Add a :class:`.CompoundModel` object or ID to this set, or add multiple at
        once when ``other`` is a :class:`.CompoundSet` or :class:`.IngredientSet`"""

        match other:
            case CompoundModel():
                return CompoundSet(
                    CompoundModel.objects.filter(
                        Q(pk__in=self._queryset) | Q(pk__in=other._queryset)
                    ),
                    sort=False,
                )

            case int():
                return CompoundSet(
                    CompoundModel.objects.filter(
                        Q(pk__in=self._queryset) | Q(pk__in=other._queryset)
                    ),
                    sort=False,
                )

            case CompoundSet():
                return CompoundSet(
                    CompoundModel.objects.filter(
                        Q(pk__in=self._queryset) | Q(pk__in=other._queryset)
                    ),
                    sort=False,
                )

            case IngredientSet():
                return CompoundSet(
                    CompoundModel.objects.filter(
                        Q(pk__in=self._queryset) | Q(pk__in=other._queryset)
                    ),
                    sort=False,
                )

            case _:
                raise NotImplementedError

    def __and__(self, other: 'CompoundSet'):
        """AND set operation, returns only compounds in both sets"""

        match other:
            case CompoundSet():
                return CompoundSet(
                    CompoundModel.objects.filter(
                        Q(pk__in=self._queryset) & Q(pk__in=other.queryset)
                    ),
                    sort=False,
                )

            case _:
                raise NotImplementedError

    def __or__(self, other: 'CompoundSet'):
        """OR set operation, returns union of both sets"""

        match other:
            case CompoundSet():
                return CompoundSet(
                    CompoundModel.objects.filter(
                        Q(pk__in=self._queryset) | Q(pk__in=other.queryset)
                    ),
                    sort=False,
                )

            case _:
                raise NotImplementedError

    def __xor__(self, other: 'CompoundSet'):
        """Exclusive OR set operation, returns all compounds in either set but
        not both"""

        match other:
            case CompoundSet():
                return CompoundSet(
                    CompoundModel.objects.filter(
                        Q(Q(pk__in=self._queryset) | Q(pk__in=other.queryset))
                        & ~Q(Q(pk__in=self._queryset) & Q(pk__in=other.queryset))
                    ),
                    sort=False,
                )

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

        if not isinstance(scaffold, int):
            assert scaffold._table == 'compound'
            scaffold = scaffold.id

        values = self.db.select_where(
            query='scaffold_superstructure',
            table='scaffold',
            key=(
                f'scaffold_base = {scaffold}'
                f' AND scaffold_superstructure IN {self.str_ids}'
            ),
            multiple=True,
            none=none,
        )
        ids = [v for (v,) in values if v]

        if not ids:
            return None
        return CompoundSet(self.db, ids)

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

        variances = self.db.execute(
            f"""
        WITH nums AS (
            SELECT scaffold_base AS base, scaffold_superstructure AS elab,
            {self.db.COMPOUND_PROPERTY_FUNCTIONS['num_heavy_atoms']}(c2.compound_mol)
            - {self.db.COMPOUND_PROPERTY_FUNCTIONS['num_heavy_atoms']}(c1.compound_mol)
            AS diff
            FROM {self.db.SQL_SCHEMA_PREFIX}scaffold
            INNER JOIN {self.db.SQL_SCHEMA_PREFIX}compound AS c1
            ON scaffold_base = c1.compound_id
            INNER JOIN {self.db.SQL_SCHEMA_PREFIX}compound AS c2
            ON scaffold_superstructure = c2.compound_id
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

        variances = [v for (v,) in variances]

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

    def summary(self, return_df: bool = False) -> None:
        """Print a summary of this compound set"""

        mrich.header(self)

        from pandas import DataFrame

        sql = f"""
        SELECT tag_name,
        COUNT(DISTINCT tag_compound)
        FROM {self.db.SQL_SCHEMA_PREFIX}tag
        WHERE tag_compound IN {self.str_ids}
        GROUP BY tag_name
        ORDER BY tag_name
        """

        cursor = self.db.execute(sql)

        data = [dict(tag=a, num_compounds=b) for a, b in cursor.fetchall()]

        df = DataFrame(data)
        df = df.set_index('tag')

        # poses

        sql = f"""
        SELECT tag_name,
        COUNT(DISTINCT tag_pose)
        FROM {self.db.SQL_SCHEMA_PREFIX}tag
        INNER JOIN {self.db.SQL_SCHEMA_PREFIX}pose
        ON pose_id = tag_pose
        WHERE pose_compound IN {self.str_ids}
        GROUP BY tag_name
        ORDER BY tag_name
        """

        cursor = self.db.execute(sql)

        for tag, count in cursor.fetchall():
            df.loc[tag, 'num_poses'] = count

        # compounds with poses

        sql = f"""
        SELECT tag_name, COUNT(DISTINCT pose_compound)
        FROM {self.db.SQL_SCHEMA_PREFIX}tag
        INNER JOIN {self.db.SQL_SCHEMA_PREFIX}pose
        ON tag_pose = pose_id
        WHERE pose_compound IN {self.str_ids}
        GROUP BY tag_name
        ORDER BY tag_name
        """

        cursor = self.db.execute(sql)

        for tag, count in cursor.fetchall():
            df.loc[tag, 'num_posed_compounds'] = count

        df.loc['TOTAL', 'num_compounds'] = len(self)
        df.loc['TOTAL', 'num_poses'] = self.num_poses
        df.loc['TOTAL', 'num_posed_compounds'] = len(self.poses.compounds)

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

        from pandas import DataFrame

        sql = f"""
        SELECT tag_name,
        COUNT(DISTINCT tag_compound)
        FROM {self.db.SQL_SCHEMA_PREFIX}tag
        WHERE tag_compound IN {self.str_ids}
        GROUP BY tag_name
        ORDER BY tag_name;
        """

        cursor = self.db.execute(sql)

        data = [dict(tag=a, num_compounds=b) for a, b in cursor.fetchall()]

        df = DataFrame(data)
        df = df.set_index('tag')

        df = df.astype(int)

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

    def shuffled(self) -> 'CompoundSet':
        """Returns a randomised copy of this set"""
        copy = self.copy()
        copy.shuffle()
        return copy

    def pop(self) -> CompoundModel:
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

        data = []

        query = ['compound_id']

        if smiles:
            query.append('compound_smiles')

        if inchikey:
            query.append('compound_inchikey')

        if alias:
            query.append('compound_alias')

        if mol:
            query.append('mol_to_binary_mol(compound_mol)')

        if metadata:
            query.append('compound_metadata')

        query = ', '.join(query)

        sql = f"""
        SELECT {query}
        FROM {self.db.SQL_SCHEMA_PREFIX}compound
        WHERE compound_id IN {self.str_ids}
        """

        if debug:
            mrich.debug('querying...')
        records = self.db.execute(sql).fetchall()

        if debug:
            generator = mrich.track(records)
        else:
            generator = records

        for row in generator:
            row = list(row)

            d = dict(id=row.pop(0))

            if smiles:
                d['smiles'] = row.pop(0)

            if inchikey:
                d['inchikey'] = row.pop(0)

            if alias:
                d['alias'] = row.pop(0)

            if mol:
                d['mol'] = Mol(row.pop(0))

            if metadata and (meta_str := row.pop(0)):
                meta_dict = loads(meta_str)

                if expand_metadata:
                    for k, v in meta_dict.items():
                        d[k] = v

                else:
                    d['metadata'] = meta_dict

            data.append(d)

        df = DataFrame(data)

        if poses or num_poses:
            if debug:
                mrich.debug('adding pose column')

            lookup = self.db.get_compound_id_pose_ids_dict(self)
            if poses:
                df['poses'] = df['id'].apply(lambda x: lookup.get(x, {}))
            if num_poses:
                df['num_poses'] = df['id'].apply(lambda x: len(lookup.get(x, {})))

        if num_reactant or num_reactions:
            if debug:
                mrich.debug('adding reaction columns')
            tuples = self.db.get_reactant_product_tuples(self.ids, deduplicated=False)

            if num_reactant:
                lookup = {}
                for r, p in tuples:
                    lookup.setdefault(r, 0)
                    lookup[r] += 1
                df['num_reactant'] = df['id'].apply(lambda x: lookup.get(x, 0))

            if num_reactions:
                lookup = {}
                for r, p in tuples:
                    lookup.setdefault(p, 0)
                    lookup[p] += 1
                df['num_reactions'] = df['id'].apply(lambda x: lookup.get(x, 0))

        if scaffolds or elabs:
            if debug:
                mrich.debug('adding scaffold columns')
            tuples = self.db.get_scaffold_tuples(self.ids)

            if scaffolds:
                lookup = {}
                for b, e in tuples:
                    lookup.setdefault(e, set())
                    lookup[e].add(b)
                df['scaffolds'] = df['id'].apply(lambda x: lookup.get(x, set()))

            if elabs:
                lookup = {}
                for b, e in tuples:
                    lookup.setdefault(b, set())
                    lookup[b].add(e)
                df['elabs'] = df['id'].apply(lambda x: lookup.get(x, set()))

        if tags:
            if debug:
                mrich.debug('adding tag column')
            lookup = self.db.get_compound_tag_dict()
            df['tags'] = df['id'].apply(lambda x: lookup.get(x, {}))

        if routes:
            if debug:
                mrich.debug('adding route column')
            lookup = self.db.get_product_id_routes_dict()
            df['routes'] = df['id'].apply(lambda x: lookup.get(x, {}))

        df = df.set_index('id')

        return df

    def get_quoted(
        self,
        *,
        supplier: str = 'any',
    ) -> 'CompoundSet':
        """Get all member compounds that have a quote from given supplier

        :param supplier: supplier name (Default value = 'any')

        """

        if supplier == 'any':
            key = f'quote_compound IN {self.str_ids}'
        else:
            key = f'quote_compound IN {self.str_ids} AND quote_supplier = "{supplier}"'

        ids = self.db.select_where(
            table='quote',
            query='DISTINCT quote_compound',
            key=key,
            multiple=True,
        )

        ids = [i for (i,) in ids]
        return CompoundSet(self.db, ids)

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
        return dict(db=str(self.db.path.resolve()), indices=self.indices)

    def write_smiles_csv(
        self, file: str, tags: bool = True, split_tags: bool = True
    ) -> None:
        """Write a CSV of the smiles contained in this set to a file

        :param file: path of the CSV file
        :param tags: include tags in output
        :param split_tags: split tags into separate columns

        """
        from pandas import DataFrame

        if tags:
            records = self.db.select_where(
                table='tag',
                query='tag_compound, tag_name',
                key=f'tag_compound IN {self.str_ids}',
                multiple=True,
                none='quiet',
            )
            TAGS = {}
            if records:
                for compound_id, tag_name in records:
                    if compound_id not in TAGS:
                        TAGS[compound_id] = set()
                    TAGS[compound_id].add(tag_name)

        records = self.db.select_where(
            table=self.table,
            query='compound_id, compound_smiles',
            key=f'compound_id IN {self.str_ids}',
            multiple=True,
        )

        data = [dict(id=id, smiles=smiles) for id, smiles in records]

        if tags:
            for d in data:
                tagset = TAGS.get(d['id'], set())

                if split_tags:
                    for tag in tagset:
                        d[tag] = True
                else:
                    d['tags'] = tagset

        df = DataFrame(data)
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
        from designdb.components.recipe import Recipe

        file = str(Path(file).resolve())

        rows = []

        for r_id in mrich.track(self.reaction_ids, prefix='Solving compound recipes'):
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
                    'target-names': str(product.compound),
                    'no-steps': 0,
                    'concentration-required-mM': None,
                    'amount-required-uL': None,
                    'batch-tag': None,
                }

                for i, reaction in enumerate(sub_recipe.reactions):
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

        for i in self.indices:
            self.db.insert_tag(name=tag, compound=i, commit=False)

        mrich.print(f'Tagged {self} w/ "{tag}"')

        self.db.commit()

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
        return IngredientSet.from_compounds(
            compounds=self, amount=amount, supplier=supplier
        )

    def split_by_scaffolds(self) -> 'dict[CompoundSet, CompoundSet]':
        """Split this set into subsets clustered by scaffold compound"""

        cluster_dict = self.db.get_compound_cluster_dict(cset=self)

        subsets = {}
        for cluster, elabs in cluster_dict.items():
            cluster = CompoundSet(self.db, list(cluster))
            subsets[cluster] = CompoundSet(self.db, list(elabs))

        return subsets

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

        for i, c in mrich.track(enumerate(self), total=len(self)):
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
        result = self.db.select_where(
            query='compound_alias',
            table='compound',
            key=f'compound_id in {self.str_ids}',
            multiple=True,
        )
        return [q for (q,) in result]

    @property
    def smiles(self) -> list[str]:
        """Returns the smiles of child compounds"""
        result = self.db.select_where(
            query='compound_smiles',
            table='compound',
            key=f'compound_id in {self.str_ids}',
            multiple=True,
        )
        return [q for (q,) in result]

    @property
    def mols(self) -> 'list[Chem.Mol]':
        """Returns the molecules of child compounds"""
        from rdkit.Chem import Mol

        result = self.db.select_where(
            query='mol_to_binary_mol(compound_mol)',
            table='compound',
            key=f'compound_id in {self.str_ids}',
            multiple=True,
        )
        return [Mol(q) for (q,) in result]

    @property
    def inchikeys(self) -> list[str]:
        """Returns the inchikeys of compounds in this set"""
        result = self.db.select_where(
            query='compound_inchikey',
            table='compound',
            key=f'compound_id in {self.str_ids}',
            multiple=True,
        )
        return [q for (q,) in result]

    @property
    def tags(self) -> set[str]:
        """Returns the set of unique tags present in this compound set"""
        values = self.db.select_where(
            table='tag',
            query='DISTINCT tag_name',
            key=f'tag_compound in {self.str_ids}',
            multiple=True,
        )
        if not values:
            return set()
        return set(v for (v,) in values)

    @property
    def num_poses(self) -> int:
        """Count the poses associated to this set of compounds"""

        return self.db.count_where(table='pose', key=f'pose_compound in {self.str_ids}')

    @property
    def poses(self) -> 'PoseSet':
        """Get the poses associated to this set of compounds"""
        from .pose import PoseSet

        ids = self.db.select_where(
            query='pose_id',
            table='pose',
            key=f'pose_compound in {self.str_ids}',
            multiple=True,
            none='warning',
        )

        if not ids:
            return PoseSet(self.db, {})

        ids = [v for (v,) in ids]
        return PoseSet(self.db, ids)

    @property
    def best_placed_poses(self) -> 'PoseSet':
        """Get the best placed pose for each compound in this set"""
        from .pose import PoseSet

        query = self.db.select_where(
            table='pose',
            query='pose_id, MIN(pose_distance_score)',
            key=f'pose_compound in {self.str_ids} GROUP BY pose_compound',
            multiple=True,
        )
        ids = [i for i, s in query]
        return PoseSet(self.db, ids)

    @property
    def str_ids(self) -> str:
        """Return an SQL formatted tuple string of the :class:`.CompoundModel` IDs"""
        return str(tuple(self.ids)).replace(',)', ')')

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
        """Get a dictionary with atomtypes as keys and corresponding
        quantities/counts as values"""
        from molparse.atomtypes import combine_atomtype_dicts

        atomtype_dicts = [c.atomtype_dict for c in self]
        return combine_atomtype_dicts(atomtype_dicts)

    @property
    def num_atoms_added(self) -> list[int]:
        """Calculate the number of atoms added w.r.t the scaffold

        :returns: list of number of atoms added values

        """

        nha = self.db.COMPOUND_PROPERTY_FUNCTIONS['num_heavy_atoms']
        sql = f"""
        WITH nums AS (
            SELECT
                A.compound_id AS comp_id,
                {nha}(A.compound_mol)
                - {nha}(B.compound_mol)
                AS diff
            FROM {self.db.SQL_SCHEMA_PREFIX}compound A,
            {self.db.SQL_SCHEMA_PREFIX}compound B
            WHERE A.compound_base = B.compound_id
            AND A.compound_id IN {self.str_ids}
        )

        SELECT compound_id, diff FROM {self.db.SQL_SCHEMA_PREFIX}compound
        LEFT JOIN nums
        ON comp_id = compound_id
        WHERE compound_id IN {self.str_ids}
        """

        query = self.db.execute(sql).fetchall()

        lookup = {k: v for k, v in query}

        return [lookup[i] for i in self.indices]

    @property
    def avg_num_atoms_added(self) -> float:
        """Calculate the average number of atoms added w.r.t the scaffold

        :returns: average number of atoms added values for compounds which have a
            scaffold

        """
        nha = self.db.COMPOUND_PROPERTY_FUNCTIONS['num_heavy_atoms']
        sql = f"""
        WITH nums AS (
            SELECT
                A.compound_id AS comp_id,
                {nha}(A.compound_mol)
                - {nha}(B.compound_mol)
                AS diff
            FROM {self.db.SQL_SCHEMA_PREFIX}compound A,
            {self.db.SQL_SCHEMA_PREFIX}compound B
            WHERE A.compound_base = B.compound_id
            AND A.compound_id IN {self.str_ids}
        )

        SELECT compound_id, diff FROM {self.db.SQL_SCHEMA_PREFIX}compound
        INNER JOIN nums
        ON comp_id = compound_id
        WHERE compound_id IN {self.str_ids}
        """

        (avg,) = self.db.execute().fetchone()

        return avg

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

        sql = f"""
        SELECT COUNT(1) FROM {self.db.SQL_SCHEMA_PREFIX}scaffold
        WHERE scaffold_superstructure IN {self.str_ids}
        GROUP BY scaffold_base
        """

        counts = self.db.execute(sql).fetchall()

        counts = [c for (c,) in counts]  # + [0 for _ in range(len(self)-len(counts))]

        from hirsch import hirsch

        return hirsch(counts)

        # return -std(counts)

    @property
    def num_scaffolds_elaborated(self) -> int:
        """Count the number of scaffold compounds that have at least one elaboration in
        this set

        :returns: number of scaffold compounds

        """

        (count,) = self.db.execute(
            f"""
                SELECT COUNT(DISTINCT scaffold_base)
                FROM {self.db.SQL_SCHEMA_PREFIX}scaffold
                WHERE scaffold_superstructure IN {self.str_ids}
            """
        ).fetchone()

        return count

    @property
    def scaffolds(self) -> 'CompoundSet':
        """Get the scaffold compounds that have at least one elaboration in this set

        :returns: :class:`.CompoundSet`

        """
        return CompoundSet(self.db, self.scaffold_ids)

    @property
    def scaffold_ids(self) -> list[int]:
        """Return a list of :class:`.CompoundModel` ID's for scaffolds of this set"""
        scaffold_ids = self.db.execute(
            f"""
                SELECT DISTINCT scaffold_base FROM {self.db.SQL_SCHEMA_PREFIX}scaffold
                WHERE scaffold_superstructure IN {self.str_ids}
            """
        ).fetchall()
        return [i for (i,) in scaffold_ids]

    @property
    def num_scaffolds(self) -> int:
        """Return a count of scaffolds of this set"""
        (count,) = self.db.execute(
            f"""
                SELECT COUNT(DISTINCT scaffold_base)
                FROM {self.db.SQL_SCHEMA_PREFIX}scaffold
                WHERE scaffold_superstructure IN {self.str_ids}
            """
        ).fetchone()
        return count

    @property
    def elabs(self) -> 'CompoundSet':
        """Returns a :class:`.CompoundSet` of all compounds that are a an elaboration
        of an existing scaffold"""

        ids = self.db.select_where(
            query='scaffold_superstructure',
            table='scaffold',
            key=(
                f'scaffold_superstructure IS NOT NULL'
                f' and scaffold_base IN {self.str_ids}'
            ),
            multiple=True,
            none='quiet',
        )

        if not ids:
            return None

        ids = [q for (q,) in ids]
        return CompoundSet(self.db, ids)

    @property
    def num_elabs(self) -> int:
        """Return a count of elaborations of this set"""
        (count,) = self.db.execute(
            f"""
                SELECT COUNT(DISTINCT scaffold_superstructure)
                FROM {self.db.SQL_SCHEMA_PREFIX}scaffold
                WHERE scaffold_base IN {self.str_ids}
            """
        ).fetchone()
        return count

    @property
    def elab_df(self) -> 'pd.DataFrame':
        """Get a DataFrame summarising the elaborations in this CompoundSet"""
        from pandas import DataFrame

        cluster_dict = self.db.get_compound_cluster_dict(max_scaffolds=1)

        data = []
        for scaffold, elabs in cluster_dict.items():
            scaffold = self.db.get_compound(id=scaffold[0])
            elabs = CompoundSet(self.db, indices=elabs)
            data.append(
                dict(
                    scaffold_id=scaffold.id,
                    scaffold_compound=scaffold,
                    elabs=elabs,
                    num_elabs=len(elabs),
                )
            )

        return DataFrame(data)

    @property
    def id_num_poses_dict(self) -> dict[int, int]:
        """Get a dictionary mapping compound ids to the number of poses"""

        sql = f"""
            SELECT pose_compound, COUNT(1) FROM {self.db.SQL_SCHEMA_PREFIX}pose
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
        """Returns a list of :class:`.ReactionModel` IDs that result in members of
        this set"""
        records = self.db.select_where(
            table='reaction',
            query='reaction_id',
            key=f'reaction_product IN {self.str_ids}',
            multiple=True,
        )
        if not records:
            return None
        return [r for (r,) in records]


class IngredientSet:
    """An :class:`.Ingredient` is a :class:`.CompoundModel` with a fixed quanitity and
    an attached quote, the :class:`.IngredientSet` is a object representing multiple
    ingredients.

    .. attention::

            :class:`.IngredientSet` objects should not be created directly. Instead they
            are returned by several methods when working with :doc:`quoting` and
            :doc:`rgen`.

    Selecting ingredients in the set
    ================================

    The :class:`.IngredientSet` can be indexed like a Python list:

    ::

            ingredient = ingredient_set[0] # first ingredient

    To get the ingredient for a specific :class:`.CompoundModel` ID:

    ::

            ingredient = ingredient_set(compound_id=13)

    """

    _columns = [
        'compound_id',
        'amount',
        'quote_id',
        'supplier',
        'max_lead_time',
        'quoted_amount',
    ]

    def __init__(
        self,
        ingredients: 'None | list[Ingredient]' = None,
        supplier: str | list | None = None,
        debug: bool = False,
    ) -> None:
        """IngredientSet initialisation"""

        ingredients = ingredients or []

        self._data = DataFrame(columns=self._columns, dtype=object)

        if debug:
            mrich.debug(self._data)

        self._supplier = supplier

        for ingredient in ingredients:
            self.add(ingredient)

        for col in self._columns:
            assert col in self._data.columns, f'{col} not in df.columns'

        if debug:
            mrich.debug(self._data)

    ### DUNDERS

    def __len__(self):
        """The number of ingredients in this set"""
        return len(self._data)

    def __str__(self) -> str:
        """Unformatted string representation"""
        return f'{{Ingredient × {len(self)}}}'

    def __repr__(self) -> str:
        """ANSI ormatted string representation"""
        return f'{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}'

    def __rich__(self) -> str:
        """Representation for mrich"""
        return f'[bold underline]{self}'

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

    def __getitem__(self, key: int) -> 'Ingredient':
        """Get a member by it's index"""
        match key:
            case int():
                series = self.df.loc[key]
                return self._get_ingredient(series)

            case _:
                raise NotImplementedError

    def __iter__(self):
        """Iterate through the ingredients"""
        return iter(self._get_ingredient(s) for i, s in self.df.iterrows())

    def __call__(
        self,
        *,
        compound_id: int | None = None,
        tag: str | None = None,
    ) -> 'IngredientSet | Ingredient | CompoundSet':
        """Get members based on a compound_id or tag"""

        if compound_id:
            # get the ingredient with the matching compound ID
            matches = self.df[self.df['compound_id'] == compound_id]

            if len(matches) == 0:
                return None

            elif len(matches) != 1:
                mrich.warning(f'Multiple ingredients in set with {compound_id=}')
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

    def __contains__(self, other: CompoundModel | Ingredient | int):
        """Check if compound or ingredient is a member of this set"""
        match other:
            case CompoundModel():
                id = other.id
            case Ingredient():
                id = other.compound_id
            case int():
                id = other

        return id in set(self.compound_ids)

    @classmethod
    def from_ingredient_df(
        cls,
        df: 'DataFrame',
        supplier: str | list | None = None,
    ) -> 'IngredientSet':
        """Create an :class:`.IngredientSet` from a DataFrame

        :param db: HIPPO Database
        :param df: DataFrame of Ingredients
        :param supplier: supplier to use for all quoting, (Default value = None)

        """
        # from numpy import nan
        self = cls.__new__(cls)

        for col in cls._columns:
            if col not in df.columns:
                raise Exception(f'{col} not in df.columns')
                df[col] = None

        self._data = df.copy()
        self._supplier = supplier

        return self

    @classmethod
    def from_json(
        cls,
        path: None | str,
        supplier: str | list | None = None,
        data: None | dict = None,
    ) -> 'IngredientSet':
        """Create an :class:`.IngredientSet` from JSON data or a JSON file

        :param db: HIPPO Database
        :param path: path to JSON data (can be ``None`` if ``data`` provided)
        :param supplier: supplier to use for all quoting, (Default value = ``None``)
        :param data: optional JSON data to parse, (Default value = ``None``)

        """

        if not data:
            data = json.load(open(path))

        df = DataFrame(columns=cls._columns, dtype=object)

        for col in cls._columns:
            df[col] = data[col]

        return cls.from_ingredient_df(df=df, supplier=supplier)

    @classmethod
    def from_ingredient_dicts(
        cls,
        dicts: list[dict],
        supplier: str | list | None = None,
    ) -> 'IngredientSet':
        """Create an :class:`.IngredientSet` from :class:`.Ingredient` dictionaries

        :param db: HIPPO Database
        :param dicts: List of individual ingredient dictionaries
        :param supplier: supplier to use for all quoting, (Default value = ``None``)

        """

        df = DataFrame(dicts, dtype=object)
        return cls.from_ingredient_df(df=df, supplier=supplier)

    @classmethod
    def from_compounds(
        cls,
        *,
        compounds: 'CompoundSet | None' = None,
        ids: list[int] | None = None,
        amount: float | list[float] = 1,
        supplier: str | list | None = None,
    ) -> 'IngredientSet':
        """Create an :class:`.IngredientSet` from a :class:`.CompoundSet` or IDs

        :param compounds: :class:`.CompoundSet` to use, if ``None`` must provide
            ``ids`` and ``db`` (Default value = None)
        :param ids: CompoundModel IDs (Default value = None)
        :param db: HIPPO Database (Default value = None)
        :param amount: Amount(s) in ``mg`` (Default value = 1)
        :param supplier: supplier to use for all quoting, (Default value = ``None``)

        """

        if not ids:
            ids = compounds.ids

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

        return cls.from_ingredient_df(df)

    ### METHODS

    def get_price(
        self, supplier: str | list[str] = None, none: str = 'error', debug: bool = False
    ) -> 'Price':
        """Calculate the price with a given supplier

        :param supplier: supplier to use for all quoting, (Default value = ``None``)

        """

        pairs = {i: q for i, q in enumerate(self.df['quote_id'])}

        quote_ids = [q for q in pairs.values() if q is not None and not isnan(q)]

        if debug:
            mrich.debug('quote_ids', quote_ids)

        if quote_ids:
            qs = CataloguePriceModel.objects.filter(pk__in=quote_ids)

            if supplier:
                qs = qs.filter(supplier=supplier)

            if qs.exists():
                prices = [
                    Price(
                        amount=k.price,
                        currency=k.currency,
                    )
                    for k in qs
                ]
                quoted = sum(prices, Price.null())
            else:
                quoted = Price.null()
                self.df['quote_id'] = None
                pairs = {i: q for i, q in enumerate(self.df['quote_id'])}

        else:
            quoted = Price.null()

        if debug:
            mrich.debug('quoted', quoted)

        unquoted = [i for i, q in pairs.items() if q is None or isnan(q)]

        unquoted_price = Price.null()

        for i in unquoted:
            ingredient = self[i]

            if debug:
                mrich.debug('unquoted', i, ingredient)

            p = ingredient.price

            unquoted_price += p

            if debug:
                mrich.debug(unquoted_price)

            quote = ingredient.quote

            if not quote:
                mrich.warning('NULL Quote:', ingredient)
                continue

            self.df.loc[i, 'quote_id'] = quote.id

            assert quote.amount

            self.df.loc[i, 'quoted_amount'] = quote.amount

        if debug:
            mrich.debug('quoted', quoted)
            mrich.debug('unquoted_price', unquoted_price)
            mrich.error('end of IngredientSet.get_price()')

        return quoted + unquoted_price

    def interactive(self, **kwargs) -> None:
        """Wrapper for :meth:`.CompoundSet.interactive`"""
        self.compounds.interactive(**kwargs)

    def add(
        self,
        ingredient: 'Ingredient | None' = None,
        *,
        compound_id: int | None = None,
        amount: float | None = None,
        quote_id: int | None = None,
        supplier: str | list[str] | None = None,
        max_lead_time: float | None = None,
        quoted_amount: float | None = None,
        debug: bool = False,
    ) -> None:
        """Add an :class:`.Ingredient` to this set

        :param ingredient: :class:`.Ingredient` to be added, if ``None`` must specify
            other parameters, (Default value = None)
        :param compound_id: :class:`.CompoundModel` ID (Default value = None)
        :param amount: amount in ``mg`` (Default value = None)
        :param quote_id: :class:`.Quote` ID (Default value = None)
        :param supplier: supplier name string or list (Default value = None)
        :param max_lead_time: maximum lead-time for quoting (in days)
            (Default value = None)
        :param quoted_amount: amount of associated :class:`.Quote`
            (Default value = None)
        :param debug: increase verbosity for debugging (Default value = False)

        """

        if ingredient:
            compound_id = ingredient.compound.pk
            amount = ingredient.amount

            q = ingredient.quote

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
            if compound_id in self._data['compound_id'].values:
                index = self._data.index[
                    self._data['compound_id'] == compound_id
                ].tolist()[0]
                self._data.loc[index, 'amount'] += amount

                # discard if the quote is no longer valid
                if (a := self.df.loc[index, 'quoted_amount']) and a < self.df.loc[
                    index, 'amount'
                ]:
                    self._data.loc[index, 'quote_id'] = None
                    self._data.loc[index, 'quoted_amount'] = None

                if debug and supplier:
                    mrich.debug('Adding to existing ingredient')
                    mrich.debug(f'{self._data.loc[index, "supplier"]=}')
                    mrich.debug(f'{supplier=}')

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
                    [self._data, addition], ignore_index=True, join='inner'
                )

                if debug:
                    mrich.out(addition)

    def _get_ingredient(
        self,
        series,
    ) -> 'Ingredient':
        """Get ingredient from one of the DataFrame rows"""

        q_id = series['quote_id']

        if isinstance(q_id, float) and isnan(q_id):
            q_id = None

        return Ingredient(
            compound=CompoundModel.objects.get(pk=series['compound_id']),
            amount=series['amount'],
            quote=q_id,
            supplier=series['supplier'],
            max_lead_time=series['max_lead_time'],
        )

    def copy(self) -> 'IngredientSet':
        """Return a copy of this :class:`.IngredientSet`"""
        return IngredientSet.from_ingredient_df(self.df, supplier=self.supplier)

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

        self.df['amount'] = amount

        # if amounts are modified the quotes should be cleared
        self.df['quote_id'] = None

        assert all(self.df['supplier'].isna()) and all(self.df['max_lead_time'].isna())

        # # update quotes
        # pairs = self.db.execute(
        #     f"""
        #     WITH matching_quotes AS (
        #         SELECT quote_id, quote_compound, MIN(quote_price)
        #         FROM {self.db.SQL_SCHEMA_PREFIX}quote
        #         WHERE quote_compound IN {self.str_compound_ids}
        #         AND quote_amount >= {amount}
        #         GROUP BY quote_compound
        #     )
        #     SELECT compound_id, quote_id FROM {self.db.SQL_SCHEMA_PREFIX}compound
        #     LEFT JOIN matching_quotes ON quote_compound = compound_id
        #     WHERE compound_id IN {self.str_compound_ids}
        # """
        # ).fetchall()

        qs = CataloguePriceModel.objects.filter(
            compound__pk__in=self.compound_ids,
            quote_amount__gte=amount,
        )

        for k in qs:
            match = self.df.index[self.df['compound_id'] == k.compound.pk][0]
            self.df.loc[match, 'quote_id'] = k.quote.pk

    def get_dict(self, data_orient: str = 'list') -> dict:
        """Get serialisable dictionary

        :param data_orient: passed to ``pandas.DataFrame.to_dict``
            (Default value = 'list')

        """
        return dict(
            supplier=self.supplier,
            data=self.df.to_dict(orient=data_orient),
        )

    def pop(self) -> Ingredient:
        """Pop the last compound in this set"""
        item = self[self.df.index[-1]]
        self.df.drop(self.df.index[-1], inplace=True)
        return item

    def shuffle(self) -> None:
        """Randomises the order of compounds in this set"""
        self._data = self.df.sample(frac=1).reset_index(drop=True)

    ### PROPERTIES

    @property
    def df(self) -> 'DataFrame':
        """Access the raw DataFrame"""
        return self._data

    @property
    def price_df(self) -> 'DataFrame':
        """DataFrame including prices"""
        df = self.df.copy()
        tuples = [(i.price, i.lead_time, i.quote.supplier) for i in self]
        df['price'] = [t[0] for t in tuples]
        df['lead_time'] = [t[1] for t in tuples]
        df['quote_supplier'] = [t[2] for t in tuples]
        return df

    @property
    def price(self) -> 'Price':
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
        self.df['supplier'] = [s] * len(self)

    @property
    def smiles(self) -> list[str]:
        """SMILES for all ingredients"""
        compound_ids = list(self.df['compound_id'])
        return CompoundModel.objects.filter(
            pk__in=compound_ids,
        ).values_list('compound_smiles', flat=True)

    @property
    def inchikeys(self) -> list[str]:
        """InChI-keys for all ingredients"""
        compound_ids = list(self.df['compound_id'])
        return CompoundModel.objects.filter(
            pk__in=compound_ids,
        ).values_list('compound_inchikeys', flat=True)

    @property
    def compound_ids(self) -> list[int]:
        """CompoundModel IDs for all ingredients"""
        return list(self.df['compound_id'].values)

    @property
    def ids(self) -> list[int]:
        """CompoundModel IDs for all ingredients"""
        return self.compound_ids

    @property
    def id_amount_pairs(self) -> list[tuple]:
        """Get a list of compound ID and amount pairs"""
        return [
            (id, amount) for id, amount in self.df[['compound_id', 'amount']].values
        ]

    @property
    def str_compound_ids(self) -> str:
        """Return an SQL formatted tuple string of the :class:`.CompoundModel` IDs"""
        return str(tuple(self.df['compound_id'].values)).replace(',)', ')')

    @property
    def compounds(self) -> 'CompoundSet':
        """:class:`.CompoundSet` of all compounds in this set"""
        return CompoundSet(self.compound_ids)

    @property
    def quote_ids(self) -> list[int]:
        """Get a list of quote ID's"""

        return [q for q in self.df['quote_id'].values if not isna(q) and q is not None]
