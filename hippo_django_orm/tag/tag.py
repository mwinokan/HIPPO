from ..models import TagModel


class Tag(TagModel):

    def __str__(self):
        return f'Tag "{self.name}"'
