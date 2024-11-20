from ..models import FeatureModel
from .feature_set import FeatureSet


class Feature(FeatureModel):

    _objects = FeatureSet.as_manager()

    def __str__(self):
        return f"{self.target} -> F{self.id}"
