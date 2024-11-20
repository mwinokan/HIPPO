from .models import ReactionModel


class Reaction(ReactionModel):

    def __str__(self):
        return f"R{self.id}"
