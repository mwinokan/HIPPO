from ..models import CompoundModel
from .compound_set import CompoundSet

# from .compound_table import CompoundTable


class Compound(CompoundModel):
    class Meta:
        app_label = "hippo"

    _objects = CompoundSet.as_manager()
