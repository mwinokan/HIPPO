from ..models import TargetModel
from .target_set import TargetSet


class Target(TargetModel):

    _objects = TargetSet.as_manager()

    ### PROPERTIES

    @property
    def features(self):
        return self._features(manager="_objects").all()

    ### DUNDERS

    def __str__(self):
        return f"T{self.id}"