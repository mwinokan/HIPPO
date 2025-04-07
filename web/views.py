from django.shortcuts import render, get_object_or_404, get_list_or_404, redirect
from django.http import HttpResponse
from django.apps import apps
from django.db.models.fields.related import ForeignKey

# from .rendertypes import ContentRenderType
from hippo.custom_models import Pose, PoseSet, PoseSetMember
import plotly.graph_objects as go
import plotly.io as pio


import mrich

### function based views


def index(request):

    # campaigns

    # from hippo.projects import Campaign, Iteration

    # # iterations = Iteration.objects.all()
    # campaigns = Campaign.objects.all()

    # iteration_content = []

    # for campaign in campaigns:

    #     targets = campaign.targets

    #     for iteration in campaign.iterations.all():
    #         iteration_content.append(
    #             dict(
    #                 name=iteration.long_name,
    #                 targets=str(targets),
    #                 status=iteration.get_status_display(),
    #                 fg_style=iteration.get_status_fg_style(),
    #                 bg_style=iteration.get_status_bg_style(),
    #             )
    #         )

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
        dict(
            model_stats=model_stats,
            # iteration_content=iteration_content
        ),
    )


def pose_sdf(request, pk: int):

    from rdkit.Chem import MolToMolBlock

    pose = Pose.objects.get(id=pk)

    if not pose:
        return HttpResponse("Unknown Pose", status=400, content_type="text/plain")

    text = MolToMolBlock(pose.mol)

    return HttpResponse(text, content_type="chemical/x-mdl-sdfile")


def pose_compare(request, pks: str):

    # Split comma-separated pks (e.g., "1,2,3") into a list
    pk_list = pks.split(",")

    # Fetch poses from the database

    poses = Pose.filter(pk__in=pk_list)

    # poses = get_list_or_404(Pose, pk__in=pk_list)

    context = {
        "poses": poses,
        "colors": ["orange", "blue", "green", "yellow", "red"],
        "num_poses": len(poses),
        "last_index": len(poses) - 1,
    }

    context["mocassin"] = poses.score_set(method="mocassin")

    return render(
        request,
        "pose_compare.html",
        context,
    )


def pose_compare_3d(request, pks: str):

    pk_list = pks.split(",")

    poses = Pose.filter(pk__in=pk_list)

    if inspirations := request.GET.get("inspirations"):
        assert len(poses) == 1

        pose = poses[0]

        poses = [i for i in pose.inspirations] + [pose]

    context = {
        "poses": poses,
    }

    return render(
        request,
        "pose_compare_3d.html",
        context,
    )


### class based views

from django.views.generic import ListView, DetailView
from django.apps import apps

### auto-generated views


class BaseDetailView(DetailView):

    _uri_param_str = "<int:pk>"

    def get_context_data(self, **kwargs):

        context = super().get_context_data(**kwargs)
        context["model_name"] = self.model._meta.model_name

        fields = []
        for field in self.model._meta.get_fields():

            field_name = field.name
            render_dict = self.object.get_field_render_type(field)
            value = None

            if not field.is_relation:
                value = getattr(self.object, field_name, None)

            elif isinstance(field, ForeignKey):

                related = getattr(self.object, field_name)

                field_name = field_name  # .removeprefix("_")  # .capitalize()
                value = related

            else:

                members = getattr(self.object, field.name).all()

                field_name = field_name  # .removeprefix("_") .capitalize()
                value = members

            if value is None:
                continue

            if render_dict.get("zero_index", False):

                if members:
                    value = members[0]
                else:
                    continue

            elif hasattr(value, "__len__") and len(value) == 0:
                continue

            fields.append((field_name, value, render_dict))

        fields = sorted(fields, key=lambda x: x[0].count("_"))
        fields = [(x[0].removeprefix("_"), x[1], x[2]) for x in fields]

        context["fields"] = fields

        return context


def generate_views_for_model(model):

    class ModelListView(ListView):
        model = model
        template_name = f"model_list.html"
        context_object_name = f"{model._meta.model_name}_list"

        def get_context_data(self, **kwargs):
            context = super().get_context_data(**kwargs)
            context["model_name"] = self.model._meta.model_name

            fields = []
            for field in self.model._meta.get_fields():

                field_name = field.name

                render_dict = self.model.get_field_render_type(self.model, field)

                if field_name in self.model._list_view_fields:
                    # if not field.is_relation and field_name in self.model._list_view_fields:
                    #     # value = getattr(self.object, field_name, None)
                    # print(field_name, render_dict)
                    fields.append((field_name, render_dict))

                else:
                    continue
                #     pass

                # elif isinstance(field, ForeignKey):

                #     related = getattr(self.object, field_name)

                #     field_name = field_name.removeprefix("_")  # .capitalize()
                #     value = related

                # else:

                #     members = getattr(self.object, field.name).all()

                #     field_name = field_name.removeprefix("_").capitalize()
                #     value = members

                # if value is None:
                #     continue

                # elif hasattr(value, "__len__") and len(value) == 0:
                #     continue

            context["fields"] = fields

            return context

    class ModelDetailView(BaseDetailView):
        model = model
        template_name = f"model_detail.html"
        context_object_name = f"{model._meta.model_name}_detail"

    return ModelListView, ModelDetailView


def get_models_for_app(app_label):
    app_config = apps.get_app_config(app_label)
    return app_config.get_models()


# Get all models in the app
app_models = get_models_for_app("hippo")

# Store generated views in a dictionary
GENERATED_VIEWS = {}
for model in app_models:

    list_view, detail_view = generate_views_for_model(model)

    GENERATED_VIEWS[model._meta.model_name] = {"list_view": list_view}

    if not model._custom_detail_view:
        GENERATED_VIEWS[model._meta.model_name]["detail_view"] = detail_view

### custom views


class PoseDetailView(BaseDetailView):
    model = Pose
    template_name = f"pose_detail.html"
    context_object_name = f"pose"

    _uri_param_str = "<str:value>"

    def get_object(self):

        value = self.kwargs.get("value")

        # Check if lookup_value is an integer (PK lookup)
        if value.isdigit():
            return get_object_or_404(Pose, pk=int(value))

        # Otherwise, assume it's a custom field (e.g., 'pose_code')
        return get_object_or_404(Pose, alias=value)


class PoseSetDetailView(BaseDetailView):
    model = PoseSet
    template_name = f"poseset_detail.html"
    context_object_name = f"poseset"

    _uri_param_str = "<int:pk>"

    def get_context_data(self, *args, **kwargs):
        context = super().get_context_data(**kwargs)

        fig = self.object.generate_tsnee_fig()
        fig = pio.to_json(fig)
        context["fig"] = fig

        return context

    def post(self, request, *args, **kwargs):
        # Get the instance you are viewing (retrieved based on pk in URL)
        obj = self.get_object()

        # Here, you can get the data from the POST request (e.g., from a form)
        new_review = request.POST.get("review", None)

        assert new_review

        pose_id = request.POST.get("pose", None)

        assert pose_id

        member = PoseSetMember.objects.get(parent_id=obj.id, pose_id=pose_id)

        assert member

        member.review = new_review

        member.save()

        # Redirect to a success URL or reload the detail view with the updated data
        return redirect("poseset_detail", pk=obj.pk)

    # def get_object(self):

    #     value = self.kwargs.get("value")

    #     # Check if lookup_value is an integer (PK lookup)
    #     if value.isdigit():
    #         return get_object_or_404(PoseSet, pk=int(value))

    #     # Otherwise, assume it's a custom field (e.g., 'pose_code')
    #     return get_object_or_404(PoseSet, alias=value)


GENERATED_VIEWS["pose"]["detail_view"] = PoseDetailView
GENERATED_VIEWS["poseset"]["detail_view"] = PoseSetDetailView
