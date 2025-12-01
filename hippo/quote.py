"""Classes to work with quote data"""

import mcol
import mrich

from .price import Price


class Quote:
    """Supplier quote for a specific quantity of a :class:`.Compound`.

    .. attention::

            :class:`.Quote` objects should not be created directly. Instead use :meth:`.Compound.get_quotes`.

    """

    _db = None

    def __init__(
        self,
        db: "Database",
        id: int,
        compound: int,
        smiles: str,
        supplier: str,
        catalogue: str,
        entry: str,
        amount: float,
        price: float,
        currency: str,
        purity: float,
        lead_time: int,
        date: str | None = None,
        type: str | None = None,
    ) -> None:
        """Quote initialisation"""

        price = Price(price, currency)

        self._db = db
        self._id = id
        self._compound = compound
        self._smiles = smiles
        self._supplier = supplier
        self._catalogue = catalogue
        self._entry = entry
        self._amount = amount
        self._price = price
        self._purity = purity
        self._lead_time = lead_time
        self._date = date
        self._type = type

        from datetime import datetime

        quote_age = (datetime.today() - datetime.strptime(self.date, "%Y-%m-%d")).days
        # if quote_age > 30:
        # mrich.warning(f'Quote is {quote_age} days old')
        # mrich.warning(self)

    ### FACTORIES

    @classmethod
    def combination(
        cls,
        required_amount: float,
        quotes: list["Quote"],
        debug: bool = False,
    ) -> "Quote":
        """Combine a list of quotes into one :class:`.Quote` object.

        * Start with biggest pack
        * Estimate by scaling linearly with unit price

        :param required_amount: amount in mg
        :param quotes: list of quotes to be combined

        """

        biggest_pack = sorted(quotes, key=lambda x: x.amount)[-1]
        unit_price = biggest_pack.price / biggest_pack.amount
        estimated_price = unit_price * required_amount

        quote_data = dict(
            db=biggest_pack.db,
            id=None,
            compound=biggest_pack.compound,
            smiles=biggest_pack.smiles,
            supplier=biggest_pack.supplier,
            catalogue=biggest_pack.catalogue,
            entry=biggest_pack.entry,
            amount=required_amount,
            price=estimated_price.amount,
            currency=biggest_pack.currency,
            purity=biggest_pack.purity,
            lead_time=biggest_pack.lead_time,
            date=biggest_pack.date,
            type=f"estimate from quote={biggest_pack.id}",
        )

        if debug:
            mrich.debug(f"Quote.combination()")
            mrich.debug(f"{required_amount=}")
            for quote in quotes:
                mrich.debug(quote)
            mrich.debug(f"{biggest_pack=}")
            mrich.debug(f"{unit_price=}")
            mrich.debug(f"{estimated_price=}")
            mrich.print(quote_data)

        self = cls.__new__(cls)
        self.__init__(**quote_data)

        return self

    ### PROPERTIES

    @property
    def entry_str(self) -> str:
        """Unformatted string including the supplier, catalogue (if available), and entry name of the quote"""
        if self.catalogue:
            return f"{self.supplier}:{self.catalogue}:{self.entry}"
        else:
            return f"{self.supplier}:{self.entry}"

    @property
    def db(self) -> "Database":
        """Returns a pointer to the parent database"""
        return self._db

    @property
    def id(self) -> int:
        """Returns the quote's database ID"""
        return self._id

    @property
    def compound(self) -> int:
        """Returns the associated :class:`.Compound`"""
        return self._compound

    @property
    def smiles(self) -> str:
        """Returns the catalogue SMILES string"""
        return self._smiles

    @property
    def supplier(self) -> str:
        """Name of the supplier"""
        return self._supplier

    @property
    def catalogue(self) -> str | None:
        """Name of the catalogue"""
        return self._catalogue

    @property
    def entry(self) -> str:
        """Name/ID of the catalogue entry"""
        return self._entry

    @property
    def amount(self) -> float:
        """Amount in mg"""
        return self._amount

    @property
    def price(self) -> "Price":
        """Price"""
        return self._price

    @property
    def currency(self) -> str:
        """Currency of the associated :class:`.Price` object"""
        return self.price.currency

    @property
    def purity(self) -> float:
        """Purity fraction"""
        return self._purity

    @property
    def lead_time(self) -> float:
        """Lead time in days"""
        return self._lead_time

    @property
    def date(self) -> str:
        """Date the quote was registered to the database"""
        return self._date

    @property
    def type(self) -> str:
        """Description of this quote"""
        return self._type

    @property
    def dict(self) -> dict:
        """Dictionary representation of this quote"""
        return dict(
            id=self.id,
            compound=self.compound,
            smiles=self.smiles,
            supplier=self.supplier,
            catalogue=self.catalogue,
            entry=self.entry,
            amount=self.amount,
            price=self.price,
            purity=self.purity,
            lead_time=self.lead_time,
            date=self.date,
            type=self.type,
        )

    @property
    def currency_symbol(self) -> str:
        """Currency symbol of the associated :class:`.Price`"""
        return self.price.symbol

    ### DUNDERS

    def __str__(self):
        """Unformatted string representation"""
        if self.purity:
            purity = f" @ {self.purity:.0%}"
        else:
            purity = ""

        if self.supplier == "Stock":
            return f"C{self.compound} In Stock: {self.amount:}mg{purity}"
        elif self.type:
            return f"C{self.compound} {self.entry_str} {self.amount:}mg{purity} = {self.price:} ({self.lead_time} days) {self.smiles} [{self.type}]"
        else:
            return f"C{self.compound} {self.entry_str} {self.amount:}mg{purity} = {self.price:} ({self.lead_time} days) {self.smiles}"

    def __repr__(self) -> str:
        """ANSI Formatted string representation"""
        return f"{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}"

    def __rich__(self) -> str:
        """Rich Formatted string representation"""
        return f"[bold underline]{self}"
