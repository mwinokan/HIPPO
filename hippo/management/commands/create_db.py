from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
import mrich


class Command(BaseCommand):
    help = "Runs makemigrations, migrate, and createsuperuser in sequence."

    def handle(self, *args, **kwargs):

        mrich.h1("create_db")

        mrich.print("Starting project setup...")

        # Step 1: Run makemigrations
        try:
            mrich.print("Running makemigrations...")
            call_command("makemigrations")
            mrich.print("Makemigrations completed.")
        except CommandError as e:
            mrich.error(f"Error during makemigrations: {e}")
            return

        # Step 2: Run migrate
        try:
            mrich.print("Running migrate...")
            call_command("migrate")
            mrich.print("Migrate completed.")
        except CommandError as e:
            mrich.error(f"Error during migrate: {e}")
            return

        # Step 3: Run createsuperuser
        try:
            mrich.print("Creating superuser...")
            call_command("createsuperuser", interactive=True)
            mrich.print("Superuser created successfully.")
        except CommandError as e:
            mrich.error(f"Error during createsuperuser: {e}")
            return

        mrich.success("Project setup complete.")
