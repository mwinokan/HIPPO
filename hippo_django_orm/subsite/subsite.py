from ..models import SubsiteModel
from .subsite_set import SubsiteSet


class Subsite(SubsiteModel):

    _objects = SubsiteSet.as_manager()

    # @property
    # def poses(self):
    #     return self._poses(manager="_objects").all()
