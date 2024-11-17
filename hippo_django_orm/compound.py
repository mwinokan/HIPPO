from .models import CompoundModel
from .compound_set import CompoundSet


class Compound(CompoundModel):

    _objects = CompoundSet.as_manager()

    def __str__(self):
        return f"C{self.id}"
