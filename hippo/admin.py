from django.contrib import admin

# import models

from .custom_models import MODELS

# setup custom wrappers

from .orm.setup import setup_models

setup_models()

# register models

for model in MODELS:
    admin.site.register(model)
