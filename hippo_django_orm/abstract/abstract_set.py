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

    ### DUNDERS

    # __getitem__ handled by models.QuerySet
