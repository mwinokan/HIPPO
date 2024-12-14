from ..models import PoseScoreModel
from .pose_score_set import PoseScoreSet


class PoseScore(PoseScoreModel):

    _objects = PoseScoreSet.as_manager()

    _module_name = "pose_score"
