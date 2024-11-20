from ..models import ReactionModel
from .reaction_set import ReactionSet


class Reaction(ReactionModel):

    _objects = ReactionSet.as_manager()

    def __str__(self):
        return f"R{self.id}"
