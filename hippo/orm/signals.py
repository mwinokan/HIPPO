from django.db.backends.signals import connection_created
from django.dispatch import receiver


@receiver(connection_created)
def load_chemicalite(connection, **kwargs) -> None:
    """Enable the chemicalite extension to sqlite3"""
    if connection.vendor != "sqlite":
        return
    connection.connection.enable_load_extension(True)
    connection.connection.load_extension("chemicalite")
