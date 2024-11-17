from .models import TargetModel


class Target(TargetModel):

    def __str__(self):
        return f"T{self.id}"
