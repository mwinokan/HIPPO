from ..models import InteractionModel
from .interaction_set import InteractionSet


class Interaction(InteractionModel):

    _objects = InteractionSet.as_manager()
