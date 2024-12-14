from ..models import CompoundScoreModel
from .compound_score_set import CompoundScoreSet


class CompoundScore(CompoundScoreModel):

    _objects = CompoundScoreSet.as_manager()
