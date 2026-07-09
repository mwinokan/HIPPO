# minimal settings to allow creating a migration for sqlite3 database,
# in case user opts to use a local sqlite database


# run like:
# cd hippo
# DJANGO_SETTINGS_MODULE=xchem_hippo.sqlite_migration_settings \
#     python -m django makemigrations designdb

from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent

INSTALLED_APPS = [
    'designdb.apps.DesigndbConfig',
]

# needs to define a database to run migrations
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        # make sure it's out of source tree
        'NAME': BASE_DIR.parent / 'dev.sqlite3',
    }
}

SECRET_KEY = 'migration-only'
DEFAULT_AUTO_FIELD = 'django.db.models.AutoField'

# Important so Django actually creates tables
MANAGE_MODELS = True
