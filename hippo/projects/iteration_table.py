from .iteration import Iteration
from ..abstract import AbstractTable


class IterationTable(AbstractTable):

    _name = "all iterations"
    _model = Iteration
