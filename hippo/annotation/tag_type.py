from ..models import TagTypeModel
from .tag_type_set import TagTypeSet


class TagType(TagTypeModel):

    _objects = TagTypeSet.as_manager()
    _module_name = "tag_type"
