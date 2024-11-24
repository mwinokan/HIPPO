import mcol
import mrich
from django.db import models
import importlib
from ..orm.managers import ManagerRouter


class AbstractModel(models.Model):
    class Meta:
        app_label = "hippo"
        abstract = True

    _id = models.BigAutoField(primary_key=True)

    def __init__(self, auto_save: bool = False, *args, **kwargs):

        #     # set up parent instance
        #     super().__setattr__("_auto_save", False)
        super().__init__(*args, **kwargs)
        self._wrapped_field_names = None

    #     super().__setattr__("_auto_save", auto_save)

    ### MODEL SETUP

    @classmethod
    def _setup_wrappers(cls):

        name = cls.__name__

        for field in dir(cls):

            mrich.debug(name, field, str(type(getattr(cls, field))))

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
    def id(self):
        return self._id

    @property
    def shorthand(self) -> str | None:
        return self._shorthand

    @property
    def wrapped_field_names(self):
        if not self._wrapped_field_names:
            self._wrapped_field_names = self.get_wrapped_field_names()
        return self._wrapped_field_names

    ### METHODS

    def get_wrapped_field_names(self) -> set[str]:

        fields = self._meta.get_fields()

        names = []
        for field in fields:
            name = field.name
            name = name.removeprefix("_")
            if name in names:
                names = [n for n in names if n != name]
            names.append(name)

        return set(names)

    def summary(self):

        from rich.panel import Panel
        from rich.box import SIMPLE_HEAVY
        from rich.table import Table

        fields = self.wrapped_field_names

        table = Table(title=self.__rich__(), box=SIMPLE_HEAVY)
        table.add_column("Field", style="var_name")
        table.add_column("Value", style="result")

        table.add_row(f"[bold]Model", f"[bold var_type]{self.__name__}")

        for field in fields:
            value = getattr(self, field)

            if not value:
                continue

            if hasattr(value, "__rich__"):
                s = value.__rich__()
            else:
                s = str(value)
            table.add_row(f"[bold]{field}", s)

        panel = Panel(table, expand=False)
        mrich.print(panel)

    # def save(self, full_clean: bool = True, *args, **kwargs):
    #     if full_clean:
    #         self.full_clean()
    #     super().save(*args, **kwargs)

    # def save(self, full_clean: bool = True, *args, **kwargs):
    #     auto = self._auto_save
    #     super().__setattr__("_auto_save", False)
    #     # self.full_clean()
    #     super().save(*args, **kwargs)
    #     if auto:
    #         super().__setattr__("_auto_save", True)

    ### DUNDERS

    def __str__(self) -> str:

        if sh := self.shorthand:
            s = f"{self.shorthand}{self.id}"
        else:
            s = f"{self.__name__}_{self.id}"

        if field := self._name_field:
            return f'{s} "{getattr(self, field)}"'

        # if hasattr(self, "name") and (name := self.name):
        # return f'{s} "{name}"'

        return s

    def __repr__(self) -> str:
        """ANSI formatted string representation"""
        return f"{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}"

    def __rich__(self) -> str:
        """Representation for mrich"""
        return f"[bold]{self}"

    @property
    def __name__(self):
        return self.__class__.__name__

    # def __setattr__(self, key, value):
    #     if self._auto_save and key in self.wrapped_field_names:
    #         result = super().__setattr__(key, value)
    #         mrich.debug(f"auto-saving on __setattr__({key=})")
    #         self.save()
    #         return result

    #     return super().__setattr__(key, value)

    # def __super_setattr__(self, key, value):
    #     return super().__setattr__(key, value)
