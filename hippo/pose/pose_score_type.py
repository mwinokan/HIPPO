from ..models import PoseScoreTypeModel
from .pose_score_type_set import PoseScoreTypeSet


class PoseScoreType(PoseScoreTypeModel):

    _objects = PoseScoreTypeSet.as_manager()
    _module_name = "pose_score_type"
    _parent_module = "pose"
