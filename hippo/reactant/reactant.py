from ..models import ReactantModel
from .reactant_set import ReactantSet


class Reactant(ReactantModel):

    _objects = ReactantSet.as_manager()

    def __str__(self):
        return f"{self.reaction} -> Reactant_{self.id}"
