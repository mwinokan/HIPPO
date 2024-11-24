import mrich
import django
from django.conf import settings

DJANGO_SETUP = False


def setup_django(
    databases,
    force: bool = False,
    include_contenttypes: bool = False,
    include_templates: bool = False,
    debug: bool = False,
):

    global DJANGO_SETUP
    if DJANGO_SETUP and not force:
        mrich.debug("django already setup")

    # define DATABASES dictionary
    DATABASES = {"default": {}}
    for alias, path in databases.items():
        DATABASES[alias] = {
            "ENGINE": "django.db.backends.sqlite3",
            "NAME": path,
            "APP_LABEL": "hippo",
        }

    INSTALLED_APPS = ["hippo"]
    TEMPLATES = []

    if include_templates:

        from pathlib import Path

        BASE_DIR = Path(__file__).resolve().parent.parent

        include_contenttypes = True
        INSTALLED_APPS.append("django.contrib.admin")
        INSTALLED_APPS.append("django.contrib.auth")
        INSTALLED_APPS.append("django.contrib.sessions")
        INSTALLED_APPS.append("django.contrib.messages")

        TEMPLATES = [
            {
                "BACKEND": "django.template.backends.django.DjangoTemplates",
                "DIRS": [
                    BASE_DIR / "templates"
                ],  # Specify the path to your template files
                "APP_DIRS": True,  # Enables template loading from `templates` directory within apps
                "OPTIONS": {
                    "context_processors": [
                        "django.template.context_processors.debug",
                        "django.template.context_processors.request",
                        "django.contrib.auth.context_processors.auth",
                        "django.contrib.messages.context_processors.messages",
                    ],
                },
            },
        ]

    if include_contenttypes:
        INSTALLED_APPS.append("django.contrib.contenttypes")

    SETTINGS = dict(
        DATABASES=DATABASES,
        INSTALLED_APPS=INSTALLED_APPS,
        TEMPLATES=TEMPLATES,
    )

    if debug:
        SETTINGS["DEBUG"] = True
        SETTINGS["LOGGING"] = {
            "version": 1,
            "filters": {
                "require_debug_true": {
                    "()": "django.utils.log.RequireDebugTrue",
                }
            },
            "handlers": {
                "console": {
                    "level": "DEBUG",
                    "filters": ["require_debug_true"],
                    "class": "logging.StreamHandler",
                }
            },
            "loggers": {
                "django.db.backends": {
                    "level": "DEBUG",
                    "handlers": ["console"],
                }
            },
        }

    # Update the settings with the custom DATABASES dictionary
    settings.configure(**SETTINGS)

    # Initialize Django
    django.setup()
    DJANGO_SETUP = True

    mrich.debug("django setup")
