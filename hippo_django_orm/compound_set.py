from django.db import models


class CompoundSet(models.QuerySet):
    ...

    def __str__(self):
        return "{" f"C x {len(self)}" "}"
