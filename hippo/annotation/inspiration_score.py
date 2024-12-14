from ..models import InspirationScoreModel
from .inspiration_score_set import InspirationScoreSet


class InspirationScore(InspirationScoreModel):

    _objects = InspirationScoreSet.as_manager()
    _module_name = "inspiration_score"
    _parent_module = "annotation"
