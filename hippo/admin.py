from django.contrib import admin

from .models import MODELS

for model in MODELS:
    admin.site.register(model)
