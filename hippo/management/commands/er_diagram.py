from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
import mrich
from django.apps import apps


class Command(BaseCommand):
    help = "Runs makemigrations, migrate, and createsuperuser in sequence."

    def handle(self, *args, **kwargs):

        mrich.h1("er_diagram")

        app_label = "hippo"
        output_file = "docs/source/hippo_er_diagram.dot"

        graph_models_args = {
            "outputfile": output_file,
            "inheritance": False,  # Adjust verbosity if needed
            "exclude_models": "*Model",  # Adjust verbosity if needed
            "verbosity": 1,  # Adjust verbosity if needed
        }

        # Step 1: Run makemigrations
        try:
            mrich.print("Running graph_models")
            mrich.writing(output_file)
            call_command("graph_models", app_label, **graph_models_args)
        except CommandError as e:
            mrich.error(f"Error during graph_models: {e}")
            return

        mrich.success("visualise on dreampuf.github.io/GraphvizOnline")
