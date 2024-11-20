from .models import PoseModel


class Pose(PoseModel):

    def __str__(self):
        return f"{self.compound} -> P{self.id}"
