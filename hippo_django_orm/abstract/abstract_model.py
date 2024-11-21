import mrich
from django.db import models


class AbstractModel(models.Model):
    class Meta:
        app_label = "hippo"
        abstract = True

    id = models.BigAutoField(primary_key=True)

    @property
    def shorthand(self) -> str | None:
        return self._shorthand

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

                    # mrich.warning(f"No property wrapper for {name}.{prop} defined")
                    # mrich.debug(f"Created property wrapper for {name}.{prop}")

                else:
                    mrich.warning(f"Existing property wrapper for {name}.{prop}")

    @classmethod
    def _create_reverse_wrapper(cls, field):

        def wrap_reverse_manager(self):
            return getattr(self, f"_{field}")(manager="_objects").all()

        setattr(cls, field, property(wrap_reverse_manager))

    def __str__(self) -> str:
        if s := self.shorthand:
            return f"{self.shorthand}{self.id}"
        else:
            return f"{self.__name__}_{self.id}"
