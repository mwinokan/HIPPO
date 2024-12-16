from ..models import IterationModel
from .iteration_set import IterationSet


class Iteration(IterationModel):

    _objects = IterationSet.as_manager()
    _parent_module = "projects"

    @property
    def long_name(self):
        return f"{self.campaign.name}.{self.number}"

    def get_status_bg_style(self):

        return {
            "P": "--color-status-red",
            "D": "--color-status-medium",
            "M": "--color-status-medium",
            "T": "--color-status-medium",
            "A": "--color-status-medium",
            "F": "--color-status-good",
            "C": "--color-status-neutral",
        }[self.status]

    def get_status_fg_style(self):

        return {
            "P": "--text-color-status-red",
            "D": "--text-color-status-medium",
            "M": "--text-color-status-medium",
            "T": "--text-color-status-medium",
            "A": "--text-color-status-medium",
            "F": "--text-color-status-good",
            "C": "--text-color-status-neutral",
        }[self.status]
