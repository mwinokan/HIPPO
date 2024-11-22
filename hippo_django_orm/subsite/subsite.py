from ..models import SubsiteModel
from .subsite_set import SubsiteSet


class Subsite(SubsiteModel):

    _objects = SubsiteSet.as_manager()

    # _poses = models.ManyToManyField(
    #     "Pose",
    #     through="ObservationModel",
    #     through_fields=("subsite", "pose"),
    #     related_name="_subsites",
    # )
