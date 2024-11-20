from ..models import TargetModel
from .target_set import TargetSet


class Target(TargetModel):

    _objects = TargetSet.as_manager()

    def __str__(self):
        return f"T{self.id}"
