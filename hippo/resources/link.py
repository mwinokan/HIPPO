from ..models import LinkModel
from .link_set import LinkSet
from pathlib import Path


class Link(LinkModel):

    _objects = LinkSet.as_manager()
    _parent_module = "resources"
