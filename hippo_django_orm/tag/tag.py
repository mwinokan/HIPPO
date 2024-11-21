from ..models import TagModel
from .tag_set import TagSet


class Tag(TagModel):

    _objects = TagSet.as_manager()

    def __str__(self):
        return f'Tag "{self.name}"'
