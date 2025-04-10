import mcol

import mrich

CURRENCIES = {
    "USD": "$",
    "EUR": "€",
    "GBP": "£",
}


class Price:
    """Class to represent a certain amount of currency. Supported currencies:

    ::

            CURRENCIES = {
                    'USD':'$',
                    'EUR':'€',
                    'GBP':'£',
            }

    """

    def __init__(self, amount, currency):

        if currency not in CURRENCIES:
            assert currency is None, f"Unrecognised {currency=}"
            assert not amount, f"Null Price can't have {amount=}"
            amount = None

        self._amount = amount
        self._currency = currency

    ### FACTORIES

    @classmethod
    def null(cls) -> "Price":
        """Zero in any currency"""
        self = cls.__new__(cls)
        self.__init__(None, None)
        return self

    @classmethod
    def from_dict(
        cls,
        d: dict,
    ) -> "Price":
        """Create a :class:`.Price` object from a dictionary:

        ::

                dict(amount: float, currency: str)

        :param d: dictionary in the above format:

        """
        self = cls.__new__(cls)
        self.__init__(d["amount"], d["currency"])
        return self

    ### PROPERTIES

    @property
    def symbol(self) -> str:
        """Currency symbol"""
        return CURRENCIES[self.currency]

    @property
    def currency(self) -> str:
        """Currency string"""
        return self._currency

    @property
    def amount(self) -> float:
        """Amount"""
        return self._amount

    @property
    def is_null(self) -> bool:
        """Is this :meth:`.Price.null` or zero?"""
        return self.amount is None

    ### METHODS

    def get_dict(self) -> dict:
        """Dictionary in the format:

        ::

                dict(amount: float, currency: str)

        """
        return dict(amount=self.amount, currency=self.currency)

    def copy(self) -> "Price":
        """Return a copy of this :class:`.Price`"""
        return Price(amount=self.amount, currency=self.currency)

    ### DUNDERS

    def __str__(self) -> str:
        """Unformatted string representation"""
        if self.currency is None:
            return "Null Price"

        return f"{self.symbol}{self.amount:.2f} {self.currency}"

    def __repr__(self) -> str:
        """ANSI Formatted string representation"""
        return f"{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}"

    def __rich__(self) -> str:
        """Rich Formatted string representation"""
        return f"[bold underline]{self}"

    def __add__(self, other: "Price") -> "Price":
        """Add two :class:`.Price` objects

        :param other: :class:`.Price` object
        :returns: :class:`.Price` object

        """

        if other is None:
            return self

        if other.is_null:
            return self

        if self.is_null:
            return other

        if self.currency != other.currency:
            raise NotImplementedError(
                f"Adding two different currencies: {self.currency} != {other.currency}"
            )
        return Price(self.amount + other.amount, self.currency)

    def __truediv__(self, other: "Price | float | int") -> "Price | float":
        """Divide this :class:`.Price` by another object

        :param other: :class:`.Price` or float or int
        :returns: :class:`.Price` object or float

        """

        if isinstance(other, int) or isinstance(other, float):
            if self.is_null:
                return self
            return Price(amount=self.amount / other, currency=self.currency)

        elif isinstance(other, Price):
            assert self.currency == other.currency
            assert not other.is_null
            return self.amount / other.amount

        raise TypeError(f"Division not supported between Price and {type(other)}")

    def __mul__(self, other: "Price | float | int") -> "Price | float":
        """Multiply this :class:`.Price` by another object

        :param other: :class:`.Price` or float or int
        :returns: :class:`.Price` object or float

        """

        if isinstance(other, int) or isinstance(other, float):
            if self.is_null:
                return self
            return Price(amount=self.amount * other, currency=self.currency)

        raise TypeError(f"Multiplication not supported between Price and {type(other)}")

    def __eq__(self, other: "Price") -> bool:
        """Compare two :class:`.Price` objects"""

        if isinstance(other, int) or isinstance(other, float):
            if self.is_null:
                return other == 0
            return self.amount == other

        if self.is_null and other.is_null:
            return True

        if self.is_null and not other.is_null:
            return False

        if not self.is_null and other.is_null:
            return False

        assert (
            self.currency == other.currency
        ), f"Comparing different currencies: {self.currency} != {other.currency}"
        return self.amount == other.amount

    def __lt__(self, other: "Price") -> bool:
        """Compare two :class:`.Price` objects"""

        if isinstance(other, int) or isinstance(other, float):
            if self.is_null:
                return False
            return self.amount > other

        if self.is_null and other.is_null:
            return False

        if self.is_null and not other.is_null:
            return True

        if not self.is_null and other.is_null:
            return False

        assert (
            self.currency == other.currency
        ), f"Comparing different currencies: {self.currency} != {other.currency}"
        return self.amount < other.amount

    def __gt__(self, other: "Price") -> bool:
        """Compare two :class:`.Price` objects"""

        if isinstance(other, int) or isinstance(other, float):
            if self.is_null:
                return False
            return self.amount < other

        if self.is_null and other.is_null:
            return False

        if self.is_null and not other.is_null:
            return False

        if not self.is_null and other.is_null:
            return True

        assert (
            self.currency == other.currency
        ), f"Comparing different currencies: {self.currency} != {other.currency}"
        return self.amount > other.amount
