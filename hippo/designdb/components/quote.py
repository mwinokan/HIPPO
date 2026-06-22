"""Component wrapping a catalogue price (quote).

A :class:`.Quote` wraps a :class:`.CataloguePriceModel` row. It can also be an
*estimated* quote (see :meth:`.Quote.estimate`), not backed by a saved row, for
when no single catalogue pack covers the required amount.
"""

import mcol
from designdb.models import CataloguePriceModel

from .price import Price


class Quote:
    """A catalogue price for a :class:`.CompoundModel`, with a fixed amount.

    Wraps a :class:`.CataloguePriceModel`. Estimated quotes (built by
    :meth:`.estimate`) wrap an *unsaved* model instance and therefore have
    ``id is None`` and :attr:`.is_estimate` ``True``.
    """

    def __init__(self, instance: CataloguePriceModel):
        """Quote initialisation"""
        self._instance = instance

    ### FACTORIES

    @classmethod
    def estimate(cls, required_amount: float, quotes: 'list[Quote]') -> 'Quote | None':
        """Estimate a quote for ``required_amount`` when no single pack is big enough.

        Scales the biggest available pack's unit price to the required amount; the
        returned quote wraps an *unsaved* :class:`.CataloguePriceModel` (``id is None``).

        :param required_amount: amount in ``mg``
        :param quotes: available :class:`.Quote` packs to scale from
        :returns: the estimated :class:`.Quote`, or ``None`` if there's nothing usable
            to scale from
        """

        usable = [q for q in quotes if q.amount and q.price is not None]

        if not usable:
            return None

        biggest_pack = max(usable, key=lambda q: q.amount)

        unit_price = biggest_pack.price / biggest_pack.amount
        estimated_price = unit_price * required_amount

        instance = CataloguePriceModel(
            catalogue_compound_id=biggest_pack._instance.catalogue_compound_id,
            vendor=biggest_pack._instance.vendor,
            supplier=biggest_pack.supplier,
            amount=required_amount,
            price=estimated_price,
            currency=biggest_pack.currency,
            purity=biggest_pack._instance.purity,
            lead_time=biggest_pack.lead_time,
        )

        return cls(instance)

    ### PROPERTIES

    @property
    def id(self) -> int | None:
        """Database ID of the wrapped quote, or ``None`` for an estimate"""
        return self._instance.pk

    @property
    def is_estimate(self) -> bool:
        """``True`` if this quote is an estimate not backed by a saved row"""
        return self._instance.pk is None

    @property
    def price(self) -> float | None:
        """Price amount (currency-less float), see :attr:`.as_price`"""
        return self._instance.price

    @property
    def currency(self) -> str | None:
        """Currency of the price"""
        return self._instance.currency

    @property
    def as_price(self) -> 'Price':
        """The price as a :class:`.Price` object"""
        return Price(self._instance.price, self._instance.currency)

    @property
    def amount(self) -> float | None:
        """Quoted amount in ``mg``"""
        return self._instance.amount

    @property
    def supplier(self) -> str | None:
        """Supplier of the quote"""
        return self._instance.supplier

    @property
    def vendor(self) -> str | None:
        """Vendor of the quote"""
        return self._instance.vendor

    @property
    def lead_time(self) -> int | None:
        """Lead time of the quote (in days)"""
        return self._instance.lead_time

    ### DUNDERS

    def __str__(self) -> str:
        """Plain string representation"""
        tag = 'estimate' if self.is_estimate else f'Q{self.id}'
        return f'{tag}: {self.as_price} for {self.amount}mg'

    def __repr__(self) -> str:
        """ANSI formatted string representation"""
        return f'{mcol.bold}{mcol.underline}{str(self)}{mcol.unbold}{mcol.ununderline}'

    def __rich__(self) -> str:
        """Representation for mrich"""
        return f'[bold underline]{str(self)}'
