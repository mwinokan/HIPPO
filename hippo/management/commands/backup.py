from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
from pathlib import Path
import mrich
import sqlite3
from ...tools import dt_hash


class Command(BaseCommand):
    help = "Backup the SQLite database to a specified directory."

    def add_arguments(self, parser):
        # Add an optional argument for the backup directory
        parser.add_argument(
            "--backup-dir",
            type=str,
            default="db/",
            help='Directory to store the database backup (default: "backups")',
        )
        parser.add_argument(
            "--pages",
            type=int,
            default=10_000,
            help="Number of pages to copy at once",
        )

    def handle(self, *args, **kwargs):

        mrich.h1("backup")

        # Get the database file location from Django settings
        db_path = Path(settings.DATABASES["default"].get("NAME"))
        if not db_path or not db_path.exists():
            raise CommandError("Could not find the SQLite database file.")

        # Get the backup directory
        backup_dir = Path(kwargs["backup_dir"])
        pages = kwargs["pages"]

        backup_dir.mkdir(exist_ok=True)

        timestamp = dt_hash()

        backup_filename = db_path.name.replace(".sqlite", f".backup.{timestamp}.sqlite")
        destination = backup_dir / backup_filename

        with mrich.spinner(f"Backing up..."):
            mrich.var("source", db_path)
            mrich.var("target", destination)

            mrich.writing(destination)

            def progress(status, remaining, total):
                mrich.debug(f"Copied {total-remaining} of {total} pages...")

            src = sqlite3.connect(db_path)
            dst = sqlite3.connect(destination)

            with dst:
                src.backup(dst, pages=pages, progress=progress)

            src.close()
            dst.close()

        mrich.success("Backed up", db_path.name)
