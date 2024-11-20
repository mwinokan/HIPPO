from ..models import PoseModel
from .pose_set import PoseSet


class Pose(PoseModel):

    _objects = PoseSet.as_manager()

    def __str__(self):
        return f"P{self.id}"
