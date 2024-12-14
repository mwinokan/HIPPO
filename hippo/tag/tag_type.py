from ..models import TagTypeModel
from .tag_type_set import TagTypeSet


class TagType(TagTypeModel):

    _objects = TagTypeSet.as_manager()

    def __str__(self):
        return f'"{self.name}"'
