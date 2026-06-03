import os
import sys
from pathlib import Path

import django
import mrich
from django.conf import settings

from .ta_auth_connector import get_auth_target_access

# fix path
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))


def configure_django(db_config, manage_models: bool):

    if settings.configured:
        return

    if manage_models:
        # sqlite3 db, create and manage models
        database = {
            'ENGINE': 'django.db.backends.sqlite3',
            'NAME': db_config,
        }
    else:
        # postgres, existing installation, don't touch

        database = {
            'ENGINE': 'django.db.backends.postgresql',
            'NAME': db_config['DB_NAME'],
            'USER': db_config['DB_USER'],
            'PASSWORD': db_config['DB_PASSWORD'],
            'HOST': db_config['DB_HOST'],
            'PORT': db_config['POSTGRES_PORT'],
            'OPTIONS': {
                # sets the schema
                'options': '-c search_path=rdkit,designdb'
            },
        }

    settings.configure(
        INSTALLED_APPS=[
            'designdb.apps.DesigndbConfig',
        ],
        DATABASES={'default': database},
        SECRET_KEY='runtime',
        DEFAULT_AUTO_FIELD='django.db.models.BigAutoField',
        TIME_ZONE='UTC',
        USE_TZ=True,
        MIGRATION_MODULES={'designdb': None},
        MANAGE_MODELS=manage_models,
    )

    django.setup()


def load_hippo(
    *,
    target_name: str,
    target_access_string: str,
    username: str,
    db: str | Path | dict | None = None,
    # copy_from: str | Path | None = None,
    # overwrite_existing: bool = False,
    # update_legacy: bool = False,
):
    """Initialisation function for HIPPO object.

    User should not call HIPPO directly because the db needs to be initialised.
    """

    mrich.bold('Creating HIPPO animal')
    mrich.var('target_name', target_name, color='arg')

    tas_list = get_auth_target_access(username)

    # mock response until auth pod is externally accessible
    tas_list = ('lb18145-1')

    if not target_access_string in tas_list:
        mrich.error(f'User {username} does not have access to {target_access_string}')
        return

    if db is None:
        # populate from env

        db = {
            'DB_NAME': os.environ.get('DB_NAME', ''),
            'DB_USER': os.environ.get('DB_USER', ''),
            'DB_PASSWORD': os.environ.get('DB_PASSWORD', ''),
            'DB_HOST': os.environ.get('DB_HOST', ''),
            'POSTGRES_PORT': os.environ.get('POSTGRES_PORT', '5432'),
        }

    if isinstance(db, str):
        # sqlite db

        db_path = Path(db)

        mrich.var('db_path', db_path, color='file')

        # if copy_from:
        #     self._db = Database.copy_from(
        #         source=copy_from,
        #         destination=db_path,
        #         animal=self,
        #         update_legacy=update_legacy,
        #         overwrite_existing=overwrite_existing,
        #     )
        # else:
        #     self._db = Database(db_path, animal=self, update_legacy=update_legacy)

        configure_django(db_path, manage_models=True)

        from django.apps import apps
        from django.db import connection

        with connection.schema_editor() as schema_editor:
            for model in apps.get_models():
                if model._meta.managed:
                    schema_editor.create_model(model)

    else:
        # postgres db
        # pass

        # self._db = PostgresDatabase(animal=self, **db)
        configure_django(db, manage_models=False)

        # import .testmodule
    from designdb.animal import HIPPO

    animal = HIPPO(target_name, target_access_string)

    mrich.success('Initialised animal', f'{target_name}')
    return animal
