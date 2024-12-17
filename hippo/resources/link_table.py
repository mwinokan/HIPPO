from .link import Link
from ..abstract import AbstractTable


class LinkTable(AbstractTable):

    _name = "all links"
    _model = Link
