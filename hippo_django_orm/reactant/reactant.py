from ..models import ReactantModel


class Reactant(ReactantModel):

    def __str__(self):
        return f"{self.reaction} -> Reactant_{self.id}"
