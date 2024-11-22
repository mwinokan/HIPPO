from ..models import SolventModel
from .solvent_set import SolventSet


class Solvent(SolventModel):

    _objects = SolventSet.as_manager()

    def __str__(self) -> str:
        return f'Solvent_{self.id} "{self.name}"'
