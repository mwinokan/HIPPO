from ..models import CompoundModel
from .compound_set import CompoundSet

# from .compound_table import CompoundTable


class Compound(CompoundModel):

    _objects = CompoundSet.as_manager()
