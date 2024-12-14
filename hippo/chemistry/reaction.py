from ..models import ReactionModel
from .reaction_set import ReactionSet


class Reaction(ReactionModel):

    _objects = ReactionSet.as_manager()
    _parent_module = "chemistry"
