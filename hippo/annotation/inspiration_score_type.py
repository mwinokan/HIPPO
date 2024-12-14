from ..models import InspirationScoreTypeModel
from .inspiration_score_type_set import InspirationScoreTypeSet


class InspirationScoreType(InspirationScoreTypeModel):

    _objects = InspirationScoreTypeSet.as_manager()
    _module_name = "inspiration_score_type"
    _parent_module = "annotation"
