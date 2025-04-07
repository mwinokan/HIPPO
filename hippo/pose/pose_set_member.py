from ..models import PoseSetMemberModel
from .pose_set_members import PoseSetMembers


class PoseSetMember(PoseSetMemberModel):

    _objects = PoseSetMembers.as_manager()
    _parent_module = "pose"

    _shorthand = "P"
