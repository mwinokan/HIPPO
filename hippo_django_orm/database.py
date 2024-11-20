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
        path: str | Path,
        db_alias: str = "default",
        # reference: str | Path | None = None,
        setup_django: bool = True,
        debug: bool = True,
        **kwargs,
    ) -> None:

        self._db_alias = db_alias

        if setup_django:
            from .setup import setup_django

            databases = {
                db_alias: path,
            }
            setup_django(databases, **kwargs)

        # create the tables
        self.create_tables()

        mrich.success(f'Connected Database "{db_alias}" @', f"[file]{path}")

    @classmethod
    def multi(cls, paths):

        databases = {
            "default": paths[0],
            # "reference":path2,
        }

        for i, path in enumerate(paths[1:]):
            databases[f"reference{i+1}"] = path

        from .setup import setup_django

        setup_django(databases)

        dbs = []
        for alias, path in databases.items():
            db = cls.__new__(cls)
            db.__init__(path=path, db_alias=alias, setup_django=False)
            dbs.append(db)

        return dbs

    @property
    def connection(self):
        return connections[self._db_alias]

    def create_tables(self, debug: bool = False) -> None:
        """For each table in the schema, rename the Model, create the table if it doesn't exist, and add any missing columns"""

        from .target import Target
        from .compound import Compound
        from .pose import Pose
        from .quote import Quote
        from .tag import Tag
        from .models import SubsiteMember
        from .subsite import Subsite
        from .interaction import Interaction
        from .feature import Feature
        from .solvent import Solvent

        # define the table names
        MODELS = [
            Target,
            Compound,
            Pose,
            Quote,
            Tag,
            Subsite,
            Interaction,
            Feature,
            SubsiteMember,
            Solvent,
        ]

        for model in MODELS:

            if debug:
                mrich.var("model", model)

            # create the table
            self.create_table(model)

            # add missing columns
            missing_columns = self.get_missing_columns(model)
            if missing_columns:
                raise ValueError("Database using legacy schema")
                # mrich.warning(
                #     f"Database using legacy schema, updating table '{table_name}'"
                # )
                # mrich.var("missing_columns", missing_columns)
                # self.update_table(model, debug=debug)

    def create_table(self, model) -> None:
        """Create a table if it doesnt exist"""
        with self.connection.schema_editor() as schema_editor:
            if model._meta.db_table not in self.connection.introspection.table_names():
                schema_editor.create_model(model)

    def get_missing_columns(self, model) -> set[str]:

        missing_columns = set()

        with self.connection.schema_editor() as schema_editor:
            # Check if the table exists
            if model._meta.db_table in self.connection.introspection.table_names():
                # Get the current columns in the table
                current_columns = [field.column for field in model._meta.fields]

                # Get the database columns
                database_columns = self.connection.introspection.get_table_description(
                    self.connection.cursor(), model._meta.db_table
                )
                database_column_names = [column.name for column in database_columns]

                # Check if each field in the model exists in the database table
                for field in model._meta.fields:
                    if field.column not in database_column_names:
                        missing_columns.add(field.column)

        return missing_columns

    def update_table(self, model, debug: bool = True) -> None:
        """Update table if you added fields (doesn't drop fields)"""
        with self.connection.schema_editor() as schema_editor:
            # Check if the table exists
            if model._meta.db_table in self.connection.introspection.table_names():
                # Get the current columns in the table
                current_columns = [field.column for field in model._meta.fields]

                # Get the database columns
                database_columns = self.connection.introspection.get_table_description(
                    self.connection.cursor(), model._meta.db_table
                )
                database_column_names = [column.name for column in database_columns]

                # Check if each field in the model exists in the database table
                for field in model._meta.fields:
                    if field.column not in database_column_names:
                        # Add the new column to the table
                        schema_editor.add_field(model, field)
                        if debug:
                            mrich.debug(f"Adding {field} to {model}")
