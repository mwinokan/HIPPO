
HIPPO "animal" object
=====================

.. py:function:: hippo.HIPPO(target_name, target_access_string, *, username, db=None, stack='production', auth_token=None)

    Factory function that configures Django, connects to the database, and returns a :class:`~designdb.animal.HIPPO` animal object bound to the specified target.

    Exposed as ``HIPPO`` at the package level: ``from hippo import HIPPO``.

    :param target_name: Name of the protein target
    :param target_access_string: Target access string for authentication
    :param username: Your FedID username (keyword-only)
    :param db: Database path (str/Path for SQLite) or connection dict for PostgreSQL. If ``None``, reads from environment variables.
    :param stack: Fragalysis stack to use (default: ``'production'``)
    :param auth_token: Optional authentication token

.. autoclass:: designdb.animal.HIPPO
    :members:
