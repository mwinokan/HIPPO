from ..models import SolventModel
from .solvent_set import SolventSet


class Solvent(SolventModel):

    _objects = SolventSet.as_manager()
    _parent_module = "chemistry"

    def __str__(self) -> str:
        return f'"{self.name}"'
