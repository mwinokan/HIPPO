from ..models import PoseModel
from .pose_set import PoseSet


class Pose(PoseModel):

    _objects = PoseSet.as_manager()
    _parent_module = "pose"
