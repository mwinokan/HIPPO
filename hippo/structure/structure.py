from ..models import StructureModel
from .structure_set import StructureSet


class Structure(StructureModel):

    _objects = StructureSet.as_manager()

    # def __str__(self) -> str:
    #     return f'"{self.name}"'
