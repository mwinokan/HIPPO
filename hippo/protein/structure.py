from ..models import StructureModel
from .structure_set import StructureSet


class Structure(StructureModel):

    _objects = StructureSet.as_manager()
    _parent_module = "protein"
