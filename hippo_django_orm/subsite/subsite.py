from ..models import SubsiteModel


class Subsite(SubsiteModel):

    def __str__(self):
        return f"Subsite_{self.id}"
