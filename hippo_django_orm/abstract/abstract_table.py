import mcol


class AbstractTable:
    """Defines methods common to all hippo QuerySet wrappers, e.g. CompoundSet"""

    _max_str_ids = 5

    ### FACTORIES

    ### PROPERTIES

    @property
    def name(self) -> str | None:
        """Returns the name of set"""
        return self._name

    @property
    def shorthand(self) -> str:
        """Returns the shorthand for elements in this set"""
        return self._model._shorthand

    @property
    def _all_objects(self):
        """access all Django Model instances"""
        return self._model._objects.all()

    ### METHODS

    def get(self, *args, **kwargs):
        return self._all_objects.get(*args, **kwargs)

    ### DUNDERS

    def __str__(self) -> str:
        """Unformatted string representation"""

        if self.name:
            s = f"{self.name}: "
        else:
            s = ""

        n = len(self)

        if n == 0:
            s += "{ empty " f"{self.__class__.__name__}" " }"

        elif n <= self._max_str_ids:
            s += "{ "
            s += ", ".join(str(member) for member in self)
            s += " }"

        else:
            s += "{" f"{self.shorthand} × {len(self)}" "}"

        return s

    def __repr__(self) -> str:
        """ANSI formatted string representation"""
        return f"{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}"

    def __rich__(self) -> str:
        """Representation for mrich"""
        return f"[bold underline]{self}"

    def __len__(self) -> int:
        """Total number of compounds"""
        return self._all_objects.count()

    def __iter__(self):
        """Iterate through members"""
        return iter(self._all_objects)

    def __getitem__(self, key):
        """Get member based on type of key:

        int: self.get(id=key)
        """

        match key:
            case int():
                return self.get(id=key)
            case _:
                raise NotImplementedError
