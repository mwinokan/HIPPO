from django.shortcuts import (
    render,
    get_object_or_404,
    get_list_or_404,
    redirect,
    reverse,
)
from django.http import HttpResponse, HttpResponseRedirect, JsonResponse

# from django.contrib import messages
from django.apps import apps
from django.db.models.fields.related import ForeignKey, ManyToManyField
from django.db.models.fields.reverse_related import (
    OneToOneRel,
    ManyToManyRel,
    ManyToOneRel,
)
from django.db.models import Count, Avg
from django.contrib.auth.decorators import login_required
from django.core.management import call_command
from django.contrib.staticfiles import finders

# from .rendertypes import ContentRenderType
from hippo.models import *
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
from .forms import *
import subprocess
from hippo.tools import sanitise_smiles, SanitisationError

import mrich
from mrich import print
import re


### function based views


def index(request):

    context = {}

    ### Targets

    context["targets"] = Target.objects.all().order_by("name")

    ### Model Stats

    model_stats = []
    for model in MODELS:
        if model._exclude_from_index:
            continue

        model_stats.append(
            dict(
                model_name=model.__name__,
                model_count=model.objects.count(),
                model=model,
            )
        )

    context["model_stats"] = model_stats

    # context["search_form"] = SearchForm()

    return render(
        request,
        "index.html",
        context,
    )


def pose_sdf(request, pk: int):

    from rdkit.Chem import MolToMolBlock

    pose = Pose.objects.get(id=pk)

    if not pose:
        return HttpResponse("Unknown Pose", status=400, content_type="text/plain")

    text = MolToMolBlock(pose.mol)

    return HttpResponse(text, content_type="chemical/x-mdl-sdfile")


def structure_pdb(request, pk: int):

    structure = Structure.objects.get(id=pk)

    if not structure:
        return HttpResponse("Unknown Structure", status=400, content_type="text/plain")

    text = structure.pdb_block

    return HttpResponse(text, content_type="chemical/x-pdb")


def pose_compare(request, pks: str):

    pk_list = pks.split(",")

    poses = Pose.objects.filter(pk__in=pk_list)

    # if inspirations := request.GET.get("inspirations"):
    #     assert len(poses) == 1

    #     pose = poses[0]

    #     poses = [i for i in pose.inspirations] + [pose]

    context = {
        "poses": poses,
        "colors": ["orange", "blue", "green", "yellow", "red"],
        "num_poses": len(poses),
        "last_index": len(poses) - 1,
    }

    context["pose_ids"] = [p.id for p in poses]

    # context["mocassin"] = poses.score_set(method="mocassin")

    return render(
        request,
        "pose_compare.html",
        context,
    )


def pose_compare_3d(request, pks: str):

    pk_list = pks.split(",")

    poses = Pose.objects.filter(pk__in=pk_list)

    if inspirations := request.GET.get("inspirations"):
        assert len(poses) == 1

        pose = poses[0]

        poses = [i for i in pose.inspirations.all()] + [pose]

    context = {
        "poses": poses,
        "colors": ["orange", "blue", "green", "yellow", "red"],
        "model_pills": request.GET.get("model_pills", False),
    }

    return render(
        request,
        "pose_compare_3d.html",
        context,
    )


@login_required
def fragalysis_download(request):

    if request.method == "POST":
        form = FragalysisDownloadForm(request.POST)

        if form.is_valid():
            download = form.save()  # Save to the database

            # Run the management command
            subprocess.Popen(
                ["python", "manage.py", "load_fragalysis", str(download.id)],
            )

            return HttpResponseRedirect(
                f"{reverse('fragalysis_download')}?id={download.id}"
            )

        else:
            return HttpResponseRedirect(reverse("fragalysis_download"))

    else:

        download_id = int(request.GET.get("id", 0))

        if not download_id:
            # use custom template
            form = FragalysisDownloadForm()
            return render(request, "fragalysis_download.html", dict(form=form))

        else:

            # use model_detail template
            download = FragalysisDownload.objects.get(id=download_id)
            context = get_base_detail_context(
                dict(object=download),
                FragalysisDownload,
                download,
            )
            return render(request, "model_detail.html", context)


@login_required
def sdf_upload(request):

    if request.method == "POST":
        form = SdfUploadForm(request.POST, request.FILES)

        if form.is_valid():

            uploaded_file = form.cleaned_data["input_file"]

            upload = form.save()  # Save to the database

            # Run the management command
            subprocess.Popen(
                ["python", "manage.py", "load_sdf", str(upload.id)],
            )

        return HttpResponseRedirect(f"{reverse('sdf_upload')}?id={upload.id}")

    else:

        upload_id = int(request.GET.get("id", 0))

        if not upload_id:
            # use custom template
            form = SdfUploadForm()
            return render(request, "sdf_upload.html", dict(form=form))

        else:

            # use model_detail template
            upload = SdfUpload.objects.get(id=upload_id)
            context = get_base_detail_context(
                dict(object=upload),
                SdfUpload,
                upload,
            )
            return render(request, "model_detail.html", context)


@login_required
def pill_demo(request):

    file_path = finders.find("css/style.css")

    with open(file_path, "rt") as f:

        searching = True
        styles = []
        for line in f:

            if searching:
                if "/* color by module */" in line:
                    searching = False
                continue

            match = re.search(r"--color-(.*): \w*\(.*\);", line)

            if not match:
                continue

            (style,) = match.groups()

            styles.append(style)

    pills = [
        f"""<div class="model-pill" style="
            background:var(--color-{style});
            color:var(--text-color-{style});
            ">{style}</div>"""
        for style in styles
    ]

    return render(request, "pill_demo.html", {"pills": pills})


def search(request):

    if request.method == "POST":
        form = SearchForm(request.POST)
        if form.is_valid():
            query = form.cleaned_data["query"]
            # mrich.debug("Search", query)
            # mrich.print(form.cleaned_data)

            if request.POST.get("target"):

                if match := re.match("^T([0-9]*)$", query):
                    (pk,) = match.groups()
                    # mrich.debug("target: pk", pk)
                    try:
                        target = Target.objects.get(id=pk)
                        return redirect("target_detail", pk=target.pk)
                    except Target.DoesNotExist:
                        pass

                try:
                    target = Target.objects.get(name=query)
                    return redirect("target_detail", pk=target.pk)
                except Target.DoesNotExist:
                    pass

            if request.POST.get("structure"):

                if match := re.match("^S([0-9]*)$", query):
                    (pk,) = match.groups()
                    # mrich.debug("structure: pk", pk)
                    try:
                        structure = Structure.objects.get(id=pk)
                        return redirect("structure_detail", pk=structure.pk)
                    except Structure.DoesNotExist:
                        pass

                try:
                    structure = Structure.objects.get(alias=query)
                    return redirect("structure_detail", pk=structure.pk)
                except Structure.DoesNotExist:
                    pass

            if request.POST.get("compound"):

                if match := re.match(r"^C([0-9]*)$", query):
                    (pk,) = match.groups()
                    # mrich.debug("compound: pk", pk)
                    try:
                        compound = Compound.objects.get(id=pk)
                        return redirect("compound_detail", pk=compound.pk)
                    except Compound.DoesNotExist:
                        pass

                # EXACT INCHIKEY LOOKUP
                if match := re.match(r"^[A-Z]{14}\-[A-Z]{10}(\-[A-Z])?$", query):
                    (inchikey,) = match.groups()
                    try:
                        compound = Compound.objects.get(inchikey=inchikey)
                        return redirect("compound_detail", pk=compound.pk)
                    except Compound.DoesNotExist:
                        pass

                # EXACT SMILES LOOKUP
                try:
                    smiles = sanitise_smiles(query)
                    compound = Compound.objects.get(smiles=smiles)
                    return redirect("compound_detail", pk=compound.pk)

                except Compound.DoesNotExist:
                    pass

                except SanitisationError:
                    pass

            if request.POST.get("pose"):

                if match := re.match("^P([0-9]*)$", query):
                    (pk,) = match.groups()
                    try:
                        pose = Pose.objects.get(id=pk)
                        return redirect("pose_detail", value=pose.pk)
                    except Pose.DoesNotExist:
                        pass

                # try:
                #     pose = Pose.objects.get(alias=query)
                #     return redirect('pose_detail', pk=pose.pk)
                # except Pose.DoesNotExist:
                #     pass

            form.add_error("query", "No matches found")

    else:
        form = SearchForm()

    return render(request, "search.html", {"form": form})


### class based views

from django.views.generic import ListView, DetailView
from django.apps import apps

### auto-generated views


def get_base_detail_context(context, model, instance):
    context["model_name"] = model._meta.model_name

    model_fields = model._meta.get_fields()

    # mrich.print(model_fields)

    fields = []
    for field in model_fields:

        field_name = field.name
        render_dict = instance.get_field_render_type(field)
        value = None

        if not field.is_relation:
            value = getattr(instance, field_name, None)
            render_dict.setdefault("order", 0)

        elif isinstance(field, ForeignKey):

            related = getattr(instance, field_name)
            render_dict.setdefault("order", 1)

            field_name = field_name  # .removeprefix("_")  # .capitalize()
            value = related

        elif isinstance(field, OneToOneRel):

            related = getattr(instance, field_name, None)
            render_dict.setdefault("order", 1)

            field_name = field_name  # .removeprefix("_")  # .capitalize()
            value = related

        elif (
            isinstance(field, ManyToManyRel)
            or isinstance(field, ManyToManyField)
            or isinstance(field, ManyToOneRel)
        ):

            members = getattr(instance, field.name).all()
            render_dict.setdefault("order", 2)

            field_name = field_name  # .removeprefix("_") .capitalize()
            value = members

        else:

            mrich.error("Unsupported field type", field_name, str(type(field)))

        if render_dict["content"] == "ContentRenderType.TEXT_DISPLAY":
            value = getattr(instance, f"get_{field_name}_display")()

        if value is None:
            continue

        if render_dict.get("zero_index", False):

            if members:
                value = members[0]
            else:
                continue

        if rel_name := render_dict.get("follow_related", False):
            try:
                value = getattr(value, rel_name)
            except AttributeError as e:
                value = [getattr(v, rel_name) for v in value]

        if split := render_dict.get("split", False):
            value = value.split(split)

        elif hasattr(value, "__len__") and len(value) == 0:
            continue

        fields.append((field_name, value, render_dict))

    fields = sorted(fields, key=lambda x: x[2].setdefault("order", 99))

    context["fields"] = fields

    return context


class BaseDetailView(DetailView):

    _uri_param_str = "<int:pk>"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context = get_base_detail_context(context, self.model, self.object)
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
                    fields.append((field_name, render_dict))

                else:
                    continue

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


# # Get all models in the app
# app_models = get_models_for_app("hippo")

# Store generated views in a dictionary
GENERATED_VIEWS = {}
for model in MODELS:

    list_view, detail_view = generate_views_for_model(model)

    if not model._custom_list_view:
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
        # if value.isdigit():
        return get_object_or_404(Pose, pk=int(value))

        # # Otherwise, assume it's a custom field (e.g., 'pose_code')
        # return get_object_or_404(Pose, alias=value)

    def get_context_data(self, *args, **kwargs):
        context = super().get_context_data(**kwargs)

        colors = ["orange", "blue", "green", "yellow", "red"]

        context["inspirations"] = [
            (inspiration, color)
            for inspiration, color in zip(self.object.inspirations.all(), colors)
        ]
        context["structure"] = self.object.placement.structure

        context["show_structure"] = self.request.GET.get("structure", "True") == "True"
        context["show_inspirations"] = (
            self.request.GET.get("inspirations", "False") == "True"
        )

        if self.request.user.is_authenticated:
            review = PoseReview.objects.filter(
                user=self.request.user, pose=self.object
            ).first()

            review_form = PoseReviewForm(instance=review)

            context["review_form"] = review_form

        if pset_id := self.request.GET.get("poseset_redirect", None):
            context["poseset_redirect"] = int(pset_id)

        return context

    def post(self, request, *args, **kwargs):

        pose = self.get_object()

        value = request.POST["review"]

        review = PoseReview.objects.filter(user=request.user, pose=pose).first()

        if not review:
            review = PoseReview(user=request.user, pose=pose)

        if value != review.review:
            review.review = value
            review.save()

        if pset_id := request.POST.get("poseset_redirect", None):
            return redirect("poseset_detail", int(pset_id))

        return redirect("pose_detail", pose.pk)


class PoseSetDetailView(BaseDetailView):
    model = PoseSet
    template_name = f"poseset_detail.html"
    context_object_name = f"poseset"

    _uri_param_str = "<int:pk>"

    def get_context_data(self, n_suggestions=3, *args, **kwargs):
        context = super().get_context_data(**kwargs)

        fig = self.object.generate_umap_fig()
        context["fig"] = pio.to_json(fig)

        x = list(
            self.object.poses.annotate(review_count=Count("pose__reviews")).values_list(
                "review_count", flat=True
            )
        )

        from collections import Counter

        review_count_distribution = Counter(x)

        values = []
        labels = []
        for count in range(max(review_count_distribution.keys()) + 1):
            if not count:
                labels.append(f"Unreviewed")
            else:
                labels.append(f"{count} review")
            values.append(review_count_distribution.get(count, 0))

        fig = go.Figure(
            go.Pie(values=values, labels=labels, textinfo="label+value+percent")
        )
        fig.update_layout(
            margin=dict(l=20, r=20, t=20, b=20),
            paper_bgcolor="rgba(0,0,0,0)",
            xaxis_title="#reviews",
            xaxis_range=[0, max(x)],
        )
        context["hist_review"] = pio.to_json(fig)

        x = list(
            self.object.poses.filter(predicted_score__isnull=False).values_list(
                "predicted_score", flat=True
            )
        )
        fig = go.Figure(go.Histogram(x=x))
        fig.update_layout(
            margin=dict(l=20, r=20, t=20, b=20),
            paper_bgcolor="rgba(0,0,0,0)",
            xaxis_title="Predicted Score",
            xaxis_range=[-1, 1],
        )
        context["hist_predicted"] = pio.to_json(fig)

        x = list(
            self.object.poses.filter(predicted_score__isnull=False).values_list(
                "prediction_uncertainty", flat=True
            )
        )
        fig = go.Figure(go.Histogram(x=x))
        fig.update_layout(
            margin=dict(l=20, r=20, t=20, b=20),
            paper_bgcolor="rgba(0,0,0,0)",
            xaxis_title="Prediction Uncertainty (STD)",
            xaxis_range=[0, 1],
        )
        context["hist_uncertainty"] = pio.to_json(fig)

        ### graphs of GPR dimensions not in UMAP

        data = list(
            self.object.poses.select_related("pose", "pose__reviews")
            # .filter(predicted_score__isnull=False)
            .annotate(review=Avg("pose__reviews__review")).values(
                "pose__binding_energy",
                "pose__inspiration_distance",
                "pose__ligand_energy",
                "predicted_score",
                "prediction_uncertainty",
                "review",
            )
        )

        y = [d["predicted_score"] for d in data if d["predicted_score"] is not None]
        y_err = [
            d["prediction_uncertainty"]
            for d in data
            if d["predicted_score"] is not None
        ]

        x = [
            d["pose__binding_energy"] for d in data if d["predicted_score"] is not None
        ]
        fig = go.Figure(
            go.Scatter(
                name="Predictions",
                x=x,
                y=y,
                mode="markers",
                error_y=dict(type="data", array=y_err, visible=True),
            )
        )

        fig.add_trace(
            go.Scatter(
                name="Reviews",
                x=x,
                y=[d["review"] for d in data if d["predicted_score"] is None],
                mode="markers",
                marker=dict(symbol="star"),
            )
        )
        fig.update_layout(
            margin=dict(l=20, r=20, t=20, b=20),
            paper_bgcolor="rgba(0,0,0,0)",
            xaxis_title="Binding Energy [kcal/mol]",
            yaxis_title="Predicted Score",
            yaxis_range=[-1, 1],
        )
        context["binding_energy"] = pio.to_json(fig)

        x = [
            d["pose__inspiration_distance"]
            for d in data
            if d["predicted_score"] is not None
        ]
        fig = go.Figure(
            go.Scatter(
                name="Predictions",
                x=x,
                y=y,
                mode="markers",
                error_y=dict(type="data", array=y_err, visible=True),
            )
        )

        fig.add_trace(
            go.Scatter(
                name="Reviews",
                x=x,
                y=[d["review"] for d in data if d["predicted_score"] is None],
                mode="markers",
                marker=dict(symbol="star"),
            )
        )
        fig.update_layout(
            margin=dict(l=20, r=20, t=20, b=20),
            paper_bgcolor="rgba(0,0,0,0)",
            xaxis_title="Inspiration Distance / RMSD [Å]",
            yaxis_title="Predicted Score",
            yaxis_range=[-1, 1],
        )
        context["inspiration_distance"] = pio.to_json(fig)

        x = [d["pose__ligand_energy"] for d in data if d["predicted_score"] is not None]
        fig = go.Figure(
            go.Scatter(
                name="Predictions",
                x=x,
                y=y,
                mode="markers",
                error_y=dict(type="data", array=y_err, visible=True),
            )
        )

        fig.add_trace(
            go.Scatter(
                name="Reviews",
                x=x,
                y=[d["review"] for d in data if d["predicted_score"] is None],
                mode="markers",
                marker=dict(symbol="star"),
            )
        )
        fig.update_layout(
            margin=dict(l=20, r=20, t=20, b=20),
            paper_bgcolor="rgba(0,0,0,0)",
            xaxis_title="Ligand Energy",
            yaxis_title="Predicted Score",
            yaxis_range=[-1, 1],
        )
        context["ligand_energy"] = pio.to_json(fig)

        # active learning suggestions

        unreviewed = self.object.unreviewed.filter(predicted_score__isnull=False)

        if unreviewed.count():
            context["uncertain_poses"] = unreviewed.select_related("pose").order_by(
                "-prediction_uncertainty"
            )[:n_suggestions]

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


GENERATED_VIEWS["pose"]["detail_view"] = PoseDetailView
GENERATED_VIEWS["poseset"]["detail_view"] = PoseSetDetailView
