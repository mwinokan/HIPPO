import django
from django.conf import settings

settings.configure(
    DATABASES={
        "default": {
            "ENGINE": "sqlite",
        }
    },
    INSTALLED_APPS=[
        "hippo_django_orm",
    ],
)
django.setup()
