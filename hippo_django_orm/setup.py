import mrich
import django
from django.conf import settings

DJANGO_SETUP = False


def setup_django(databases, force: bool = False):

    global DJANGO_SETUP
    if DJANGO_SETUP and not force:
        mrich.debug("django already setup")

    # define DATABASES dictionary
    DATABASES = {}
    for alias, path in databases.items():
        DATABASES[alias] = {
            "ENGINE": "django.db.backends.sqlite3",
            "NAME": path,
            "APP_LABEL": "hippo",
        }

    # Update the settings with the custom DATABASES dictionary
    settings.configure(DATABASES=DATABASES)

    # Initialize Django
    django.setup()
    DJANGO_SETUP = True

    mrich.debug("django setup")
