from django.shortcuts import render
from django.http import HttpResponse
from django.apps import apps
from django.db.models.fields.related import ForeignKey

import mrich

### function based views


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

    return render(
        request,
        "index.html",
        dict(model_stats=model_stats, iteration_content=iteration_content),
    )


### class based views

from django.views.generic import ListView, DetailView
from hippo.custom_models import Target
from django.apps import apps


def generate_views_for_model(model):

    class ModelListView(ListView):
        model = model
        template_name = f"model_list.html"
        context_object_name = f"{model._meta.model_name}_list"

        def get_context_data(self, **kwargs):
            context = super().get_context_data(**kwargs)
            context["model_name"] = self.model._meta.model_name
            return context

    class ModelDetailView(DetailView):
        model = model
        template_name = f"model_detail.html"
        context_object_name = f"{model._meta.model_name}_detail"

        def get_context_data(self, **kwargs):

            context = super().get_context_data(**kwargs)
            context["model_name"] = self.model._meta.model_name
            # context['fields'] = [field.name for field in self.model._meta.fields]

            # if hasattr(self.model, "_detail_view_skip_fields"):
            #     skip_fields = self.model._detail_view_skip_fields
            # else:
            #     skip_fields = []

            # fields = [
            #     (field.name, getattr(self.object, field.name, None))
            #     for field in self.model._meta.fields
            #     if field.name not in skip_fields and not isinstance(field, ForeignKey)
            # ]
            # context["fields"] = fields

            # fields = [
            #     (field.name, getattr(self.object, field.name, None), str(self.object.get_field_render_type(field.name)))
            #     for field in self.model._meta.fields
            #     # if field.name not in skip_fields# and not isinstance(field, ForeignKey)
            # ]
            # context["fields"] = fields

            ### related fields
            fields = []
            for field in self.model._meta.get_fields():

                field_name = field.name
                render_dict = self.object.get_field_render_type(field)
                value = None

                if not field.is_relation:
                    value = getattr(self.object, field_name, None)

                elif isinstance(field, ForeignKey):

                    related = getattr(self.object, field_name)

                    # if not related:
                    # continue

                    field_name = field_name.removeprefix("_")  # .capitalize()
                    value = related

                    mrich.print(field_name, value)

                else:

                    members = getattr(self.object, field.name).all()

                    # if not members:
                    #     continue

                    # print(members)

                    field_name = field_name.removeprefix("_").capitalize()
                    value = members

                if value is None:
                    continue

                elif hasattr(value, "__len__") and len(value) == 0:
                    continue

                fields.append((field_name, value, render_dict))

            # mrich.print(fields)

            # context["foreign_fields"] = foreign_fields
            # context["related_fields"] = related_fields
            context["fields"] = fields

            return context

    return ModelListView, ModelDetailView


# Get all models in the app
app_models = apps.get_models()

# Store generated views in a dictionary
GENERATED_VIEWS = {}
for model in app_models:
    list_view, detail_view = generate_views_for_model(model)
    GENERATED_VIEWS[model._meta.model_name] = {
        "list_view": list_view,
        "detail_view": detail_view,
    }
