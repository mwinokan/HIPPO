from django.shortcuts import render
from django.http import HttpResponse
from django.apps import apps


# Create your views here.


def index(request):

    # campaigns

    from hippo.projects import Campaign, Iteration

    # iterations = Iteration.objects.all()
    campaigns = Campaign.objects.all()

    iteration_content = []

    for campaign in campaigns:

        targets = campaign.targets

        for iteration in campaign.iterations.all():
            iteration_content.append(
                dict(
                    name=iteration.long_name,
                    targets=str(targets),
                    status=iteration.get_status_display(),
                    fg_style=iteration.get_status_fg_style(),
                    bg_style=iteration.get_status_bg_style(),
                )
            )

    # display statistics on models
    app_config = apps.get_app_config("hippo")
    models = app_config.get_models()

    model_stats = []
    for model in sorted(models, key=lambda x: x.__name__):
        try:
            model_stats.append(
                dict(model_name=model.__name__, model_count=model.objects.count())
            )
        except Exception as e:
            print(e)
            continue

    # return HttpResponse(model_stats)

    return render(
        request,
        "index.html",
        dict(model_stats=model_stats, iteration_content=iteration_content),
    )

    # return HttpResponse("Hello, world. You're at the hippo_django_app index.")
