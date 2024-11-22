import mcol
import mrich
from django.db import models
import importlib
from ..orm.managers import ManagerRouter


class AbstractModel(models.Model):
    class Meta:
        app_label = "hippo"
        abstract = True

    id = models.BigAutoField(primary_key=True)

    ### MODEL SETUP

    @classmethod
    def _check_model(cls):

        name = cls.__name__

        for field in dir(cls):

            if field == "objects":
                mrich.warning(f"No custom Manager defined for {name}")

            elif isinstance(
                getattr(cls, field),
                models.fields.related_descriptors.ReverseManyToOneDescriptor,
            ):

                if not field.startswith("_"):
                    mrich.warning(field, "ReverseManyToOneDescriptor")

                prop = field.removeprefix("_")

                if prop not in dir(cls):
                    cls._create_reverse_wrapper(prop)
                else:
                    mrich.warning(f"Existing property wrapper for {name}.{prop}")

        prop = "_table"

        if prop not in dir(cls):
            cls._create_table_wrapper()
        else:
            mrich.warning(f"Existing table wrapper: {name}.{prop}")

    @classmethod
    def _create_reverse_wrapper(cls, field) -> None:
        """Assign a custom ManagerRouter to the <field> to choose between the custom <model>Set managor of the default ManyToOne ReverseField / RelatedManager"""

        def get_manager_router(self) -> "ManagerRouter":
            return ManagerRouter(model=self, field=field)

        # Set the _<field> attribute
        setattr(cls, field, property(get_manager_router))

    @classmethod
    def _create_table_wrapper(cls) -> None:
        """Set the _table property on the child <model> Model to point to the <model>Table class if defined"""

        # define names
        module_name = cls.__name__.lower()
        model_name = cls.__name__
        table_class = f"{model_name}Table"

        # import submodule
        module = importlib.import_module(f"..{module_name}", package=__package__)

        def get_table_wrapper():
            """Get an instance of the <model>Table"""
            return getattr(module, table_class)()

        # check if the <model>Table class exists
        if hasattr(module, table_class):

            # set the _table attribute on the class
            setattr(cls, "_table", get_table_wrapper)

        else:
            mrich.warning(f"{module_name} has no {table_class}")

    ### CLASSMETHODS

    @classmethod
    def all(cls):
        if hasattr(cls, "_table"):
            # defined programmatically by _create_table_wrapper
            return cls._table()
        return cls._objects.all()

    @classmethod
    def filter(cls):
        return cls._objects.filter()

    ### PROPERTIES

    @property
    def shorthand(self) -> str | None:
        return self._shorthand

    ### METHODS

    def get_wrapped_field_names(self) -> list[str]:

        fields = self._meta.get_fields()

        names = []
        for field in fields:
            name = field.name
            name = name.removeprefix("_")
            if name in names:
                names = [n for n in names if n != name]
            names.append(name)

        return names

    def summary(self):

        from rich.table import Table

        fields = self.get_wrapped_field_names()

        table = Table(title=self.__rich__())
        table.add_column("field")
        table.add_column("value")

        for field in fields:
            value = getattr(self, field)

            if hasattr(value, "__rich__"):
                s = value.__rich__()
            else:
                s = str(value)
            table.add_row(field, s)

        mrich.print(table)

    ### DUNDERS

    def __str__(self) -> str:

        if sh := self.shorthand:
            s = f"{self.shorthand}{self.id}"
        else:
            s = f"{self.__name__}_{self.id}"

        if hasattr(self, "alias") and (alias := self.alias):
            return f'{s} "{alias}"'

        if hasattr(self, "name") and (name := self.name):
            return f'{s} "{name}"'

        return s

    def __repr__(self) -> str:
        """ANSI formatted string representation"""
        return f"{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}"

    def __rich__(self) -> str:
        """Representation for mrich"""
        return f"[bold underline]{self}"

    def __name__(self):
        return self.__class__.__name__
