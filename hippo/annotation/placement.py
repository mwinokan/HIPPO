from ..models import PlacementModel
from .placement_set import PlacementSet


class Placement(PlacementModel):

    _objects = PlacementSet.as_manager()
    _parent_module = "annotation"

    # def __str__(self) -> str:
    #     return f'"{self.name}"'
