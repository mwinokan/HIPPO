from ..models import SolventModel
from .solvent_set import SolventSet


class Solvent(SolventModel):

    _objects = SolventSet.as_manager()
