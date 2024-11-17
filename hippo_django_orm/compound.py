from .models import CompoundModel


class Compound(CompoundModel):

    def __str__(self):
        return f"C{self.id}"
