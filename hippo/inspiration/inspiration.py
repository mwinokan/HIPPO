from ..models import InspirationModel
from .inspiration_set import InspirationSet


class Inspiration(InspirationModel):

    _objects = InspirationSet.as_manager()
