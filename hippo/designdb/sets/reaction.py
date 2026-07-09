"""Classes for working with sets of :class:`.ReactionModel` objects"""

import mcol
import mrich
import pandas as pd
from designdb.models import CompoundModel, ReactantModel, ReactionModel
from designdb.sets.compound import CompoundSet
from django.db.models import Q
from IPython.display import display
from ipywidgets import (
    BoundedIntText,
    Checkbox,
    GridBox,
    Layout,
    VBox,
    interactive_output,
)


class ReactionSet:
    """Object representing a subset of the 'reaction' table in the :class:`.Database`.

    .. attention::

            :class:`.ReactionSet` objects should not be created directly. Instead use
            the :meth:`.HIPPO.reactions` property. See :doc:`getting_started` and
            :doc:`insert_elaborations`.

    Use as an iterable
    ==================

    Iterate through :class:`.ReactionModel` objects in the set:

    ::

            rset = animal.reactions[:100]

            for reaction in rset:
                    ...

    Check membership
    ================

    To determine if a :class:`.ReactionModel` is present in the set:

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

    def __init__(
        self,
        queryset=None,
        *,
        sort: bool = True,
        name: str | None = None,
    ) -> None:
        """ReactionSet initialisation"""

        if queryset:
            if isinstance(queryset, list):
                self._queryset = ReactionModel.objects.filter(pk__in=queryset)
            else:
                self._queryset = queryset
        else:
            self._queryset = ReactionModel.objects.none()

        self._name = name
        if sort:
            self._queryset = self._queryset.order_by('pk')

    def __str__(self) -> str:
        """Unformatted string representation"""

        if self.name:
            s = f'{self.name}: '
        else:
            s = ''

        s += f'{{R × {len(self)}}}'

        return s

    def __repr__(self) -> str:
        """ANSI Formatted string representation"""
        return f'{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}'

    def __rich__(self) -> str:
        """Rich Formatted string representation"""
        return f'[bold underline]{self}'

    def __len__(self) -> int:
        """Number of member :class:`.ReactionModel` objects"""
        return self._queryset.count()

    def __iter__(self):
        """Iterate through member :class:`.ReactionModel` objects"""
        return iter(self._queryset)

    def __getitem__(self, key) -> 'ReactionModel | ReactionSet':
        """Get member :class:`.ReactionModel` object by single, slice or
        list/set/tuple of ID"""

        match key:
            case int():
                try:
                    # reaction = ReactionModel.objects.get(pk=key)
                    reaction = self._queryset[key]
                except ReactionModel.DoesNotExist as exc:
                    mrich.error(f'list index out of range: {key=} for {self}')
                    raise ReactionModel.DoesNotExist from exc

                return reaction

            case slice():
                # positional slice of the (ordered) members. Slice `.all()` (a
                # fresh, unevaluated clone) so this returns a queryset even when
                # self._queryset is already evaluated (which would otherwise slice
                # to a list of instances). sort=False: can't re-order a sliced qs.
                return ReactionSet(self._queryset.all()[key], sort=False)

            case _:
                mrich.error(
                    f'Unsupported type for ReactionSet.__getitem__():'
                    f' {key=} {type(key)}'
                )

        return None

    def __add__(self, other: 'ReactionSet') -> 'ReactionSet':
        """Add a :class:`.ReactionSet` to this one"""
        if not other:
            return self.copy()
        # materialise pks so repeated accumulation (e.g. combining many recipes)
        # stays a single flat query instead of nesting pk__in subqueries, which
        # recurses and blows the recursion limit at scale
        ids = set(self._queryset.values_list('pk', flat=True)) | set(
            other.queryset.values_list('pk', flat=True)
        )
        return ReactionSet(ReactionModel.objects.filter(pk__in=ids), sort=False)

    def __sub__(
        self,
        other: 'ReactionSet',
    ) -> 'ReactionSet':
        """Substract a :class:`.ReactionSet` from this set"""
        match other:
            case ReactionSet():
                ids = set(self._queryset.values_list('pk', flat=True)) - set(
                    other.queryset.values_list('pk', flat=True)
                )
                return ReactionSet(
                    ReactionModel.objects.filter(pk__in=ids),
                    sort=False,
                )

    ### METHODS

    @classmethod
    def union(cls, sets: 'list[ReactionSet]') -> 'ReactionSet':
        """Union several :class:`.ReactionSet`\\ s into one flat set in a single pass.

        Collects pks from each set and builds one ``pk__in`` query, avoiding the
        nested subqueries that repeated ``+`` would accumulate.
        """
        ids: set[int] = set()
        for s in sets:
            ids |= set(s._queryset.values_list('pk', flat=True))
        return cls(ReactionModel.objects.filter(pk__in=ids), sort=False)

    def add(self, r: ReactionModel) -> None:
        """Add a :class:`.ReactionModel` to this set

        :param r: :class:`.ReactionModel` to be added

        """
        assert isinstance(r, ReactionModel)
        self._queryset = ReactionModel.objects.filter(
            pk__in=list(self._queryset.values_list('pk', flat=True)) + [r.pk],
        )

    def interactive(self):
        """Creates a ipywidget to interactively navigate this PoseSet."""

        a = BoundedIntText(
            value=0,
            min=0,
            max=len(self) - 1,
            step=1,
            description=f'Rs (/{len(self)}):',
            disabled=False,
        )

        b = Checkbox(description='Name', value=True)
        c = Checkbox(description='Summary', value=False)
        d = Checkbox(description='Draw', value=True)
        e = Checkbox(description='Check chemistry', value=False)
        f = Checkbox(description='ReactantModel Quotes', value=False)

        ui1 = GridBox(
            [b, c, d], layout=Layout(grid_template_columns='repeat(5, 100px)')
        )
        ui2 = GridBox([e, f], layout=Layout(grid_template_columns='repeat(2, 150px)'))
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
                'i': a,
                'name': b,
                'summary': c,
                'draw': d,
                'check_chemistry': e,
                'reactants': f,
            },
        )

        display(ui, out)

    def get_df(self, smiles=True, mols=True, **kwargs) -> pd.DataFrame:
        """Construct a pandas.DataFrame of this ReactionSet

        :param smiles: Include smiles column (Default value = True)
        :param mols: Include `rdkit.Chem.Mol` column (Default value = True)
        :param kwargs: keyword arguments are passed on to
            :meth:`.ReactionModel.get_dict:`

        """

        mrich.debug('Using slower ReactionModel.dict rather than direct SQL query...')

        data = []
        for r in mrich.track(self, prefix='ReactionSet --> DataFrame'):
            data.append(r.get_dict(smiles=smiles, mols=mols, **kwargs))

        return pd.DataFrame(data)

    def copy(self) -> 'ReactionSet':
        """Return a copy of this set"""
        return ReactionSet(self._queryset.all(), sort=False, name=self.name)

    def reverse(self) -> None:
        """Reverse the ordering of this set in-place"""
        self._queryset = self._queryset.reverse()

    def get_recipes(self, amounts: float | list[float] = 1.0, **kwargs):
        """Get the :class:`.Recipe` object(s) from this set of recipes

        :param amounts: float or list/generator of product amounts in mg,
            (Default value = 1.0)
        :param kwargs: keyword arguments are passed on to
            :meth:`.RecipeService.from_reactions`

        """
        # convenience bridge to the service layer
        from designdb.services.recipe import RecipeService

        return RecipeService.from_reactions(reactions=self, amount=amounts, **kwargs)

    def summary(self) -> None:
        """Print a summary of the Reactions"""

        mrich.header(self)
        for reaction in self:
            print(repr(reaction))

    ### PROPERTIES

    @property
    def name(self) -> str | None:
        """Returns the name of set"""
        return self._name

    @property
    def queryset(self):
        """Returns the underlying Django queryset"""
        return self._queryset

    @property
    def indices(self) -> list[int]:
        """Returns the ids of reactions in this set"""
        return list(self._queryset.values_list('pk', flat=True))

    @property
    def ids(self) -> list[int]:
        """Returns the ids of reactions in this set"""
        return list(self.indices)

    @property
    def types(self) -> list[str]:
        """Returns the unique reaction types in this set"""
        return list(self._queryset.values_list('reaction_type', flat=True).distinct())

    @property
    def num_types(self) -> int:
        """Returns the number of reaction types in this set"""
        return self._queryset.values('reaction_type').distinct().count()

    @property
    def products(self) -> CompoundSet:
        """Get all product compounds that can be synthesised with these reactions
        (no intermediates)"""

        qs = CompoundModel.objects.filter(
            pk__in=self._queryset.values('product_compound'),
        ).exclude(
            pk__in=self.intermediates.queryset.values('pk'),
        )
        cset = CompoundSet(qs)
        if self.name:
            cset._name = f'products of {self}'
        return cset

    @property
    def intermediates(self) -> CompoundSet:
        """Get all intermediate compounds that can be synthesised with these
        reactions"""

        # NB! not 100% sure about this queryset
        qs = CompoundModel.objects.filter(
            Q(
                pk__in=ReactantModel.objects.values('compound'),
            )
            & Q(pk__in=self._queryset.values('product_compound')),
        )
        cset = CompoundSet(qs)

        if self.name:
            cset._name = f'intermediates of {self}'
        return cset

    @property
    def reactants(self) -> 'CompoundSet':
        """Get all reactant compounds that are used by these reactions"""

        qs = ReactantModel.objects.filter(
            reaction__in=self._queryset,
        ).values('compound')
        cset = CompoundSet(qs)
        if self.name:
            cset._name = f'reactants of {self}'
        return cset

    @property
    def get_dict(self) -> dict[str]:
        """Serializable dictionary"""
        return dict(indices=self.indices)
