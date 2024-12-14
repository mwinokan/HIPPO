import os
import glob
from django.core.management.base import BaseCommand
from django.db import connection
from django.conf import settings
import mrich


class Command(BaseCommand):
    help = "Empty the database and delete all migration files"

    def handle(self, *args, **kwargs):

        mrich.h1("delete_db")

        mrich.h3("Starting database reset...")

        # Step 1: Empty the database
        self._delete_database()

        # Step 2: Delete migration files
        self._delete_migration_files()

        mrich.success('Database reset complete. Run "create_db" to reinitialize.')

    def _delete_database(self):
        """Delete the database"""

        db_path = settings.DATABASES["default"]["NAME"]
        mrich.warning("Deleting", db_path)
        os.remove(db_path)

    def _delete_migration_files(self):
        """Delete migration files from all apps."""
        migrations_deleted = 0
        for app in settings.INSTALLED_APPS:
            # Skip non-app entries (e.g., middleware or custom settings)
            if not os.path.exists(os.path.join(app, "migrations")):
                continue

            migrations_path = os.path.join(app, "migrations")
            for migration_file in glob.glob(
                os.path.join(migrations_path, "[0-9]*_*.py")
            ):
                mrich.warning("Deleting", migration_file)
                os.remove(migration_file)
                migrations_deleted += 1
