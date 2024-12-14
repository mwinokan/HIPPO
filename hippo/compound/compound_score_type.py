from ..models import CompoundScoreTypeModel
from .compound_score_type_set import CompoundScoreTypeSet


class CompoundScoreType(CompoundScoreTypeModel):

    _objects = CompoundScoreTypeSet.as_manager()
    _module_name = "compound_score_type"
    _parent_module = "compound"
