from ..models import SubsiteModel
from .subsite_set import SubsiteSet


class Subsite(SubsiteModel):

    _objects = SubsiteSet.as_manager()
    _parent_module = "annotation"
