from django.db import models
from .abstract_table import AbstractTable


class AbstractQuerySet(models.QuerySet, AbstractTable):

    def __init__(self, *args, **kwargs):
        models.QuerySet.__init__(self, *args, **kwargs)
        self._name = None

    ### PROPERTIES

    ### METHODS

    ### DUNDERS

    # __getitem__ handled by models.QuerySet
