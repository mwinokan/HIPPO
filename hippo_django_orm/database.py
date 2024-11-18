"""
Standalone wrapper for Django's ORM without needing any project setup or configuration

See: https://gist.github.com/mw3i/b879895272a28d1c789f23ee91555620
"""

import django
from django.conf import settings
from django.db import models, connections

import mrich
from pathlib import Path

# enable custom signal handlers
from . import signals


class Database:
    """Standalone wrapper for Django's ORM"""

    def __init__(
        self,
        # path: str | Path,
        name: str,
        debug: bool = True,
    ) -> None:

        # define DATABASES dictionary
        DATABASES = {
            "default": {
                "ENGINE": "django.db.backends.sqlite3",
                "NAME": name,
                "APP_LABEL": "hippo",
            }
        }

        # Update the settings with the custom DATABASES dictionary
        settings.configure(DATABASES=DATABASES)

        # Initialize Django
        django.setup()

        # create the tables
        self.create_tables()

        mrich.success("Connected Database @", f"[file]{name}")

    def create_tables(self, debug: bool = False) -> None:
        """For each table in the schema, rename the Model, create the table if it doesn't exist, and add any missing columns"""

        from .target import Target
        from .compound import Compound
        from .pose import Pose

        # define the table names
        MODELS = {
            "target": Target,
            "compound": Compound,
            "pose": Pose,
        }

        for table_name, model in MODELS.items():

            if debug:
                mrich.var("table_name", table_name)
                mrich.var("model", model)

            # rename the Model class
            model._meta.db_table = table_name

            # create the table
            self.create_table(model)

            # add missing columns
            missing_columns = self.get_missing_columns(model)
            if missing_columns:
                mrich.warning(
                    f"Database using legacy schema, updating table '{table_name}'"
                )
                mrich.var("missing_columns", missing_columns)
                self.update_table(model, debug=debug)

    def create_table(self, model) -> None:
        """Create a table if it doesnt exist"""
        with connections["default"].schema_editor() as schema_editor:
            if (
                model._meta.db_table
                not in connections["default"].introspection.table_names()
            ):
                schema_editor.create_model(model)

    def get_missing_columns(self, model) -> set[str]:

        missing_columns = set()

        with connections["default"].schema_editor() as schema_editor:
            # Check if the table exists
            if (
                model._meta.db_table
                in connections["default"].introspection.table_names()
            ):
                # Get the current columns in the table
                current_columns = [field.column for field in model._meta.fields]

                # Get the database columns
                database_columns = connections[
                    "default"
                ].introspection.get_table_description(
                    connections["default"].cursor(), model._meta.db_table
                )
                database_column_names = [column.name for column in database_columns]

                # Check if each field in the model exists in the database table
                for field in model._meta.fields:
                    if field.column not in database_column_names:
                        missing_columns.add(field.column)

        return missing_columns

    def update_table(self, model, debug: bool = True) -> None:
        """Update table if you added fields (doesn't drop fields)"""
        with connections["default"].schema_editor() as schema_editor:
            # Check if the table exists
            if (
                model._meta.db_table
                in connections["default"].introspection.table_names()
            ):
                # Get the current columns in the table
                current_columns = [field.column for field in model._meta.fields]

                # Get the database columns
                database_columns = connections[
                    "default"
                ].introspection.get_table_description(
                    connections["default"].cursor(), model._meta.db_table
                )
                database_column_names = [column.name for column in database_columns]

                # Check if each field in the model exists in the database table
                for field in model._meta.fields:
                    if field.column not in database_column_names:
                        # Add the new column to the table
                        schema_editor.add_field(model, field)
                        if debug:
                            mrich.debug(f"Adding {field} to {model}")
