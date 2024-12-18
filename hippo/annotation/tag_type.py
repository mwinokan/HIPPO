from ..models import TagTypeModel
from .tag_type_set import TagTypeSet


class TagType(TagTypeModel):

    _objects = TagTypeSet.as_manager()
    _parent_module = "annotation"
    _module_name = "tag_type"
