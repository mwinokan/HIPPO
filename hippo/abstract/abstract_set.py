from django.db import models
from .abstract_table import AbstractTable


class AbstractQuerySet(models.QuerySet, AbstractTable):

    _max_str_ids = 5

    def __init__(self, *args, **kwargs):
        models.QuerySet.__init__(self, *args, **kwargs)
        self._name = None

    ### PROPERTIES

    @property
    def shorthand(self):
        return self.model._shorthand

    ### METHODS

    # def filter(self, *args, **kwargs):
    #     # raise ValueError(f"{args=} {kwargs=}")
    #     # mrich.debug(self.__class__.__name__, "filter", args, kwargs)
    #     kwargs = {(f"_{k}" if not k.startswith("_") else k):v for k,v in kwargs.items()}
    #     return super().filter(*args, **kwargs)

    def filter(self, *args, **kwargs):
        # kwargs = self._add_prefix_in_kwargs(kwargs)
        return super().filter(*args, **kwargs)

    ### DUNDERS

    # __getitem__ handled by models.QuerySet
