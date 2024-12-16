import mrich

from django.db.backends.signals import connection_created

# from django.db.models.signals import pre_save
from django.dispatch import receiver


@receiver(connection_created)
def load_chemicalite(connection, **kwargs) -> None:
    """Enable the chemicalite extension to sqlite3"""
    if connection.vendor != "sqlite":
        return
    connection.connection.enable_load_extension(True)
    connection.connection.load_extension("chemicalite")
    # mrich.debug("(hippo.signals.load_chemicalite) Enabled chemicalite extension")


# @receiver(pre_save)
# def validate_pre_save(sender, instance, **kwargs):
#     instance.full_clean()  # Validate the instance
