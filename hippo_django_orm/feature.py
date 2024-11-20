from .models import FeatureModel


class Feature(FeatureModel):

    def __str__(self):
        return f"{self.target} -> F{self.id}"
