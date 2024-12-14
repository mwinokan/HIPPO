from ..models import CompoundScoreModel
from .compound_score_set import CompoundScoreSet


class CompoundScore(CompoundScoreModel):

    _objects = CompoundScoreSet.as_manager()
    _module_name = "compound_score"
    _parent_module = "compound"
