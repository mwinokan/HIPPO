from ..models import CompoundScoreTypeModel
from .compound_score_type_set import CompoundScoreTypeSet


class CompoundScoreType(CompoundScoreTypeModel):

    _objects = CompoundScoreTypeSet.as_manager()
