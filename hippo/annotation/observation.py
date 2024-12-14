from ..models import ObservationModel
from .observation_set import ObservationSet


class Observation(ObservationModel):

    _objects = ObservationSet.as_manager()
    _parent_module = "annotation"
