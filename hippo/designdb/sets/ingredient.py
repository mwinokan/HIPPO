"""IngredientSet: a set of :class:`.Ingredient`\\ s (compounds with amounts/quotes).

Builds on :class:`.CompoundSet` (the dependency runs ``ingredient -> compound``;
``CompoundSet``'s few uses of ``IngredientSet`` are deferred local imports).
"""

import json

import mcol
import mrich
from designdb.components.compound import Ingredient
from designdb.components.price import Price
from designdb.models import CataloguePriceModel, CompoundModel
from designdb.sets.compound import CompoundSet
from pandas import DataFrame, concat, isna, to_numeric


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

        for _, row in other._data.iterrows():
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
    def sum_sets(
        cls,
        sets: 'list[IngredientSet]',
        supplier: str | list | None = None,
    ) -> 'IngredientSet':
        """Merge several :class:`.IngredientSet`\\ s into one in a single pass.

        Equivalent to accumulating them with ``+=`` (amounts for a shared compound
        are summed, the first-seen quote is kept and dropped once the summed amount
        exceeds its quoted amount) but O(total ingredients) rather than the O(n^2)
        of repeated pairwise addition. See :meth:`.add`.
        """
        frames = [s._data for s in sets if not s._data.empty]
        if not frames:
            return cls(supplier=supplier)

        combined = concat(frames, ignore_index=True, join='inner')

        # first-seen row per compound keeps its quote/supplier/lead_time (even a
        # null quote) -- add() never overwrites an existing ingredient's quote
        result = combined.drop_duplicates('compound_id', keep='first').set_index(
            'compound_id'
        )

        grouped = combined.groupby('compound_id', sort=False)
        result['amount'] = grouped['amount'].sum()
        counts = grouped.size()

        # a quote is dropped only for compounds that were actually merged
        # (appear >1) whose summed amount exceeds the quoted amount (see add())
        quoted = to_numeric(result['quoted_amount'], errors='coerce')
        amount = to_numeric(result['amount'], errors='coerce')
        invalid = (
            (counts.reindex(result.index) > 1)
            & quoted.notna()
            & (quoted != 0)
            & (quoted < amount)
        )
        result.loc[invalid, 'quote_id'] = None
        result.loc[invalid, 'quoted_amount'] = None

        return cls.from_ingredient_df(result.reset_index(), supplier=supplier)

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

        # coerce to int: df values may be stored as float/object (pandas)
        quote_ids = [int(q) for q in pairs.values() if q is not None and not isna(q)]

        if debug:
            mrich.debug('quote_ids', quote_ids)

        if quote_ids:
            qs = CataloguePriceModel.objects.filter(pk__in=set(quote_ids))

            if supplier:
                qs = qs.filter(supplier=supplier)

            if qs.exists():
                # map pk -> Price, then sum over quote_ids so that ingredients
                # sharing the same catalogue row are counted with multiplicity
                # (filter(pk__in=...) collapses duplicates to one row each)
                price_by_id = {
                    k.pk: Price(amount=k.price, currency=k.currency) for k in qs
                }
                quoted = Price.null()
                for q in quote_ids:
                    price = price_by_id.get(q)
                    if price is not None:
                        quoted += price
            else:
                quoted = Price.null()
                self.df['quote_id'] = None
                pairs = {i: q for i, q in enumerate(self.df['quote_id'])}

        else:
            quoted = Price.null()

        if debug:
            mrich.debug('quoted', quoted)

        unquoted = [i for i, q in pairs.items() if q is None or isna(q)]

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

        if isinstance(q_id, float) and isna(q_id):
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
