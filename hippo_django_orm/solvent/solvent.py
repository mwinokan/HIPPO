from ..models import SolventModel


class Solvent(SolventModel):

    def __str__(self):
        return f"Solvent_{self.id}"
