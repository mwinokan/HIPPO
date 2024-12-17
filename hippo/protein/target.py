from ..models import TargetModel
from .target_set import TargetSet


class Target(TargetModel):

    _objects = TargetSet.as_manager()
    _parent_module = "protein"
