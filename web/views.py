from django.shortcuts import render
from django.http import HttpResponse
from django.apps import apps


# Create your views here.


def index(request):

    # display statistics on models
    app_config = apps.get_app_config("hippo")
    models = app_config.get_models()

    model_stats = []
    for model in models:
        try:
            model_stats.append(
                dict(model_name=model.__name__, model_count=model.objects.count())
            )
        except Exception as e:
            print(e)
            continue

    # return HttpResponse(model_stats)

    return render(request, "index.html", dict(model_stats=model_stats))

    # return HttpResponse("Hello, world. You're at the hippo_django_app index.")
