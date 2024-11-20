from .models import InteractionModel


class Interaction(InteractionModel):

    def __str__(self):
        return f"Interaction_{self.id}"
