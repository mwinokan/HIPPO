"""This file defines the base models for the HIPPO Database schema using django, these classes are then subclassed in other files to create Python wrappers"""

from pathlib import Path
from datetime import date

from django.db import models
from django.core.validators import MinValueValidator, MaxValueValidator

from .orm.expressions import MolFromSmiles, MolPatternBfpFromSmiles
from .orm.fields import MolField
from .orm.validators import validate_list_of_integers, validate_coord

from molparse.rdkit.features import FEATURE_FAMILIES, INTERACTION_TYPES
from django.shortcuts import reverse

import sys
import mrich
import mcol
from mrich import print

sys.path.append("..")
from web.rendertypes import FieldRenderType, ContentRenderType, DEFAULTS

CURRENCIES = {
    "EUR": "€",
    "USD": "$",
    "GBP": "£",
}

### ABSTRACT MODEL


class AbstractModel(models.Model):
    class Meta:
        app_label = "hippo"
        abstract = True

    id = models.BigAutoField(primary_key=True)

    _shorthand = None
    _name_field = None
    _style = "other"

    _custom_list_view = False
    _custom_detail_view = False
    _exclude_from_index = False
    _has_list_view = True

    _field_render_types = {
        "smiles": dict(
            type=FieldRenderType.TABLE,
            content=ContentRenderType.TEXT_MONOSPACE,
            copyable=True,
        ),
        "alias": dict(
            type=FieldRenderType.TABLE,
            content=ContentRenderType.TEXT_MONOSPACE,
            copyable=True,
        ),
        "origin": dict(
            type=FieldRenderType.TABLE,
            content=ContentRenderType.TEXT_MONOSPACE,
            # copyable=False,
        ),
        "name": dict(
            type=FieldRenderType.TABLE,
            content=ContentRenderType.TEXT_MONOSPACE,
            copyable=True,
        ),
        "inchikey": dict(
            type=FieldRenderType.TABLE,
            content=ContentRenderType.TEXT_MONOSPACE,
            copyable=True,
        ),
        "mol": dict(type=FieldRenderType.TABLE, content=ContentRenderType.MOL_2D_SVG),
        "pdb_block": dict(
            type=FieldRenderType.TOGGLE_CARD,
            content=ContentRenderType.TEXT_MONOSPACE,
            split="\n",
        ),
        "metadata": dict(
            type=FieldRenderType.TOGGLE_CARD, content=ContentRenderType.DICT_TABLE
        ),
    }

    _list_view_fields = ["smiles", "inchikey", "mol"]

    @property
    def model_pill_html(self):
        return f"""<div class="model-pill" style="
        background:var(--color-{self.__name__.lower()}, 
                   var(--color-{self._style}));
        color:var(--text-color-{self.__name__.lower()}, 
              var(--text-color-{self._style}));
        ">{self}</div>"""

    @property
    def class_pill_html(cls):
        return f"""<div class="model-pill" style="
        background:var(--color-{cls.__name__.lower()}, 
                   var(--color-{cls._style}));
        color:var(--text-color-{cls.__name__.lower()}, 
              var(--text-color-{cls._style}));
        ">{cls.__name__}</div>"""

    @property
    def has_list_view(self):
        return self._has_list_view

    def get_absolute_url(self):
        return reverse(f"{self.__name__.lower()}_detail", args=[str(self.id)])

    def summary(self):

        from rich.panel import Panel
        from rich.box import SIMPLE_HEAVY
        from rich.table import Table

        fields = self._meta.get_fields()

        table = Table(title=self.__rich__(), box=SIMPLE_HEAVY)
        table.add_column("Field", style="var_name")
        table.add_column("Value", style="result")

        table.add_row(f"[bold]Model", f"[bold var_type]{self.__name__}")

        for field in fields:
            value = getattr(self, field.name, None)

            if not value:
                continue

            if hasattr(value, "__rich__"):
                s = value.__rich__()
            else:
                s = str(value)
            table.add_row(f"[bold]{field.name}", s)

        panel = Panel(table, expand=False)
        mrich.print(panel)

    def clean_and_save(self):
        self.full_clean()
        self.save()

    def get_field_render_type(self, field, debug: bool = False):

        if field.name in self._field_render_types:
            data = self._field_render_types[field.name]

            if debug:
                mrich.debug("custom", field.name, type(field), data)

        else:
            data = DEFAULTS.get(str(type(field)), None)

            if debug:
                mrich.debug("default", field.name, type(field), data)

        if not data:
            data = {}

        if "type" not in data:
            mrich.warning("No default field render type", field.name, str(type(field)))
            data["type"] = "FieldRenderType.TABLE"

        if "content" not in data:
            mrich.warning(
                "No default content render type", field.name, str(type(field))
            )
            data["content"] = "ContentRenderType.TEXT_NORMAL"

        data["type"] = str(data["type"])
        data["content"] = str(data["content"])

        return data

    def get_admin_url(self):
        return reverse(f"admin:hippo_{self.__name__.lower()}_change", args=[self.id])

    ### DUNDERS

    def __str__(self) -> str:

        id_str = self.id or "?"

        if sh := self._shorthand:
            s = f"{self._shorthand}{id_str}"
        else:
            s = f"{self.__name__}_{id_str}"

        if field := self._name_field:
            name_str = getattr(self, field)
            if name_str:
                return f'{s} "{name_str}"'

        return s

    def __repr__(self) -> str:
        """ANSI formatted string representation"""
        return f"{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}"

    def __rich__(self) -> str:
        """Representation for mrich"""
        return f"[bold]{self}"

    @property
    def __name__(self):
        return self.__class__.__name__

    @property
    def class_name(self):
        return self.__class__.__name__


### PROTEIN


class Target(AbstractModel):

    name = models.CharField(max_length=200, blank=True, unique=True)

    _shorthand = "T"
    _name_field = "name"
    _style = "projects"


class Structure(AbstractModel):

    # class Meta:
    #     unique_together = ("protein_file", "alias")

    alias = models.CharField(max_length=60, blank=True, unique=True, null=True)
    pdb_block = models.TextField(blank=True, null=True)

    # files = models.ManyToManyField("File", related_name="structures")

    protein_file = models.OneToOneField(
        "File", related_name="structure", on_delete=models.PROTECT
    )

    target = models.ForeignKey(
        "Target", on_delete=models.CASCADE, related_name="structures"
    )

    STRUCTURE_ORIGINS = {
        "EXPERIMENT": "Experimentally observed structure",
        "COMPUTED": "Virtual / computed structure",
        "MANUAL": "Human-placed structure",
    }

    origin = models.CharField(
        max_length=10, choices=[(k, v) for k, v in STRUCTURE_ORIGINS.items()]
    )

    resolution = models.FloatField(blank=True, null=True)

    # poses = models.ManyToManyField(
    #     "Pose",
    #     through="Placement",
    #     through_fields=("structure", "pose"),
    #     related_name="structures",
    # )

    metadata = models.JSONField(default=dict, blank=True)

    _shorthand = "S"
    _name_field = "alias"
    _style = "protein"

    _list_view_fields = ["smiles", "target", "origin"]


### LIGAND (2D)


class Compound(AbstractModel):

    smiles = models.CharField(max_length=300, unique=True)
    inchikey = models.CharField(max_length=27, unique=True)

    # alias = models.CharField(max_length=60, blank=True, unique=True, null=True)
    # metadata = models.JSONField(default=dict, blank=True)

    mol = models.GeneratedField(
        expression=MolFromSmiles("smiles", "mol"),
        output_field=MolField(blank=True),
        db_persist=True,
    )

    pattern_bfp = models.GeneratedField(
        expression=MolPatternBfpFromSmiles("smiles", "pattern_bfp"),
        output_field=models.BinaryField(blank=True),
        db_persist=True,
    )

    # scaffolds = models.ManyToManyField(
    #     "Compound", related_name="elaborations", blank=True
    # )

    _shorthand = "C"
    # _name_field = "inchikey"
    _style = "compound"

    _field_render_types = AbstractModel._field_render_types.copy()
    _field_render_types.update(
        {
            "pattern_bfp": dict(type=FieldRenderType.HIDDEN),
            # "mol": dict(type=FieldRenderType.HIDDEN),
            # "origin": dict(
            #     type=FieldRenderType.TABLE,
            #     content=ContentRenderType.TEXT_MONOSPACE,
            #     copyable=False,
            # ),
            "tags": dict(
                type=FieldRenderType.HIDDEN,
                content=ContentRenderType.INSTANCE_PILL,
                follow_related="tag",
            ),
        }
    )

    _list_view_fields = ["smiles", "mol"]

    def get_mol_svg_text(self, width=300, height=200):

        import re
        from rdkit.Chem.Draw import MolDraw2DSVG

        drawer = MolDraw2DSVG(width, height)
        drawer.DrawMolecule(self.mol)
        drawer.FinishDrawing()
        value = drawer.GetDrawingText()

        # transparent background
        value = re.sub(r"<rect style='opacity:1.0;fill:#FFFFFF.*> <\/rect>", "", value)

        return value

    def clean(self):
        """Ensure the SMILES field is sanitized before validation."""
        if self.smiles:
            from .tools import sanitise_smiles

            self.smiles = sanitise_smiles(self.smiles)
        super().clean()

    @classmethod
    def bulk_register(
        cls,
        smiles: list[str],
        metadata: str | dict | None = None,
        radical: str = "remove",
    ) -> dict[str, str]:

        # from .compound import Compound
        from .tools import inchikey_from_smiles, sanitise_smiles
        from .orm.formatters import dict_formatter

        if metadata is None:
            metadata = [None] * len(smiles)

        # get or create instance
        objects = []
        smiles_map = {}
        for s, m in zip(smiles, metadata):

            s_new = sanitise_smiles(s, verbosity=False, radical=radical)
            i = inchikey_from_smiles(s_new)
            m = dict_formatter(m)

            smiles_map[s] = s_new

            instance = cls(
                smiles=s_new,
                inchikey=i,
            )

            objects.append(instance)

        cls.objects.bulk_create(
            objects,
            ignore_conflicts=True,
            unique_fields=["smiles", "inchikey"],
        )

        return smiles_map


### LIGAND (3D)


class Pose(AbstractModel):
    # class Meta:
    # unique_together = ("mol")

    inchikey = models.CharField(max_length=27, blank=True)
    # alias = models.CharField(max_length=60, blank=True, unique=True, null=True)
    smiles = models.CharField(max_length=300, blank=True, null=True)

    mol = MolField(blank=False, unique=True)
    metadata = models.JSONField(default=dict, blank=True)
    is_fingerprinted = models.BooleanField(default=False)

    POSE_ORIGINS = {
        "EXPERIMENT": "Experimentally observed pose",
        "COMPUTED": "Virtual / computed pose",
        "MANUAL": "Human-placed pose",
    }

    origin = models.CharField(
        max_length=10, choices=[(k, v) for k, v in POSE_ORIGINS.items()]
    )

    inspirations = models.ManyToManyField(
        "Pose",
        through="Inspiration",
        through_fields=("derivative", "original"),
        related_name="derivatives",
    )

    _shorthand = "P"
    # _name_field = "alias"
    _custom_detail_view = True

    _field_render_types = AbstractModel._field_render_types.copy()
    _field_render_types.update(
        {
            "origin": dict(
                type=FieldRenderType.TABLE,
                content=ContentRenderType.TEXT_MONOSPACE,
                copyable=False,
            ),
            "tags": dict(
                type=FieldRenderType.HIDDEN,
                content=ContentRenderType.INSTANCE_PILL,
                follow_related="tag",
            ),
        }
    )

    _list_view_fields = ["smiles"]

    def get_mol_svg_text(self, width=300, height=200):

        import re
        from rdkit.Chem.Draw import MolDraw2DSVG
        from rdkit.Chem import MolToSmiles, MolFromSmiles

        mol = MolFromSmiles(MolToSmiles(self.mol))

        drawer = MolDraw2DSVG(width, height)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        value = drawer.GetDrawingText()

        # transparent background
        value = re.sub(r"<rect style='opacity:1.0;fill:#FFFFFF.*> <\/rect>", "", value)
        value = re.sub(r"<\?xml version='1\.0' encoding='iso-8859-1'\?>", "", value)

        return value


class PoseSet(AbstractModel):

    name = models.CharField(max_length=300, unique=True)

    _name_field = "name"
    _custom_detail_view = True
    _style = "pose"

    def compute_umap(self):

        mrich.bold(self, ".compute_umap()")

        n = self.poses.count()

        mrich.var("#poses", n)

        if not n:
            mrich.warning("No poses")
            return False

        with mrich.loading("imports..."):
            from math import comb
            from itertools import combinations
            from mucos import MuCOS_score
            import numpy as np
            import umap

        distance_matrix = np.zeros((n, n))

        members = list(self.poses.prefetch_related("pose").all())

        n_combos = comb(n, 2)

        mrich.var("#pairs", n_combos)

        combos = combinations(range(n), 2)

        for i, j in mrich.track(
            combos, total=n_combos, prefix="computing distance_matrix"
        ):
            m1 = members[i].pose.mol
            m2 = members[j].pose.mol
            dist = 1 - MuCOS_score(m1, m2)
            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist

        with mrich.loading("Creating reducer"):
            reducer = umap.UMAP(metric="precomputed")

        with mrich.loading("Creating embedding"):
            embedding = reducer.fit_transform(distance_matrix)

        x = embedding[:, 0]
        y = embedding[:, 1]

        for m, x, y in zip(members, x, y):
            m.pc1 = x
            m.pc2 = y

        PoseSetMember.objects.bulk_update(members, fields=["pc1", "pc2"])

        return True

    def generate_tsnee_fig(self):

        import plotly.graph_objects as go
        from pandas import DataFrame
        import numpy as np

        # Compute UMAP embedding if needed
        if self.poses.filter(pc1__isnull=False).count() != self.poses.count():
            self.compute_umap()
        else:
            mrich.debug("No UMAP embedding needed")

        # Extrapolate human reviews
        if (
            self.poses.filter(
                review__isnull=False, predicted_score__isnull=False
            ).count()
            != 0
        ):
            self.extrapolate_reviews()
        else:
            mrich.debug("No extrapolation needed")

        members = self.poses

        df = DataFrame(
            members.values(
                "pose_id",
                "pc1",
                "pc2",
                "review",
                "predicted_score",
                "prediction_confidence",
            )
        )

        fig = go.Figure()

        if members.count():

            df_good = df[df["review"] == 1]
            df_bad = df[df["review"] == 0]
            df_test = df[df["review"].isnull()]

            if not len(df_good) and not len(df_bad):
                # no reviews, so no extrapolation
                customdata = [(pose_id, None, None) for pose_id in df["pose_id"]]
                fig.add_trace(
                    go.Scatter(
                        name="unreviewed",
                        x=df["pc1"],
                        y=df["pc2"],
                        mode="markers",
                        marker=dict(color="gray", size=10),
                        hovertemplate="%{text}<br>PC1: %{x}<br>PC2: %{y}<extra></extra>",
                        text=df["pose_id"].apply(lambda x: f"P{x}"),
                        customdata=customdata,
                    )
                )

            else:
                # GOOD HUMAN REVIEWS

                customdata = [(pose_id, 1, 1) for pose_id in df_good["pose_id"]]

                fig.add_trace(
                    go.Scatter(
                        name="GOOD",
                        x=df_good["pc1"],
                        y=df_good["pc2"],
                        mode="markers",
                        marker=dict(color="rgb(0,255,0)", symbol="star", size=10),
                        hovertemplate="%{text}<br>PC1: %{x}<br>PC2: %{y}<extra></extra>",
                        text=df_good["pose_id"].apply(lambda x: f"P{x} [GOOD]"),
                        customdata=customdata,
                    )
                )

                # BAD HUMAN REVIEWS

                customdata = [(pose_id, 0, 1) for pose_id in df_bad["pose_id"]]

                fig.add_trace(
                    go.Scatter(
                        name="BAD",
                        x=df_bad["pc1"],
                        y=df_bad["pc2"],
                        mode="markers",
                        marker=dict(color="red", symbol="star", size=10),
                        hovertemplate="%{text}<br>PC1: %{x}<br>PC2: %{y}<extra></extra>",
                        text=df_bad["pose_id"].apply(lambda x: f"P{x} [BAD]"),
                        customdata=customdata,
                    )
                )

                # EXTRAPOLATED

                data = df_test[
                    ["pose_id", "predicted_score", "prediction_confidence"]
                ].values

                customdata = [(a, b, c) for a, b, c in data]
                text = [
                    f"P{a} [GOOD?]" if b > 0.5 else f"P{a} [BAD?]" for a, b, _ in data
                ]

                fig.add_trace(
                    go.Scatter(
                        name="Predicted",
                        x=df_test["pc1"],
                        y=df_test["pc2"],
                        mode="markers",
                        marker=dict(
                            color=df_test["predicted_score"],
                            colorscale=[
                                [0.0, "red"],
                                [0.5, "white"],
                                [1.0, "rgb(0,255,0)"],
                            ],
                            opacity=np.maximum(df_test["prediction_confidence"], 0.2),
                            line=dict(color="black", width=1),
                        ),
                        hovertemplate="%{text}<br>PC1: %{x}<br>PC2: %{y}<br>Predicted: %{customdata[1]}<br>Confidence: %{customdata[2]}<extra></extra>",
                        text=text,
                        customdata=customdata,
                    )
                )

        fig.update_layout(
            margin=dict(l=20, r=20, t=20, b=20),
            paper_bgcolor="rgba(0,0,0,0)",
        )

        return fig

    def extrapolate_reviews(self):

        self.summary()

        if not self.poses.count():
            return

        from pandas import DataFrame, isna

        df = DataFrame(self.poses.values("pose_id", "pc1", "pc2", "review"))

        df_train = df[df["review"].notnull()]
        df_good = df[df["review"] == 1]
        df_bad = df[df["review"] == 0]
        df_test = df[df["review"].isnull()]

        mrich.var("#total", len(df))
        mrich.var("#train", len(df_train))
        mrich.var("#test", len(df_test))

        if not len(df_train):
            mrich.warning("Can't extrapolate, no reviews")
            return

        if not len(df_test):
            mrich.warning("Can't extrapolate, no unreviewed")
            return

        X_train = df_train[["pc1", "pc2"]].values
        y_train = df_train["review"].apply(lambda x: 1 if x else 0).values
        X_test = df_test[["pc1", "pc2"]].values

        with mrich.loading("imports..."):
            from sklearn.gaussian_process import GaussianProcessRegressor
            from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C

        kernel = C(1.0, (1e-3, 1e3)) * RBF(
            length_scale=1.0, length_scale_bounds=(1e-2, 1e2)
        )
        gpr = GaussianProcessRegressor(
            kernel=kernel, n_restarts_optimizer=10, alpha=1e-2
        )
        gpr.fit(X_train, y_train)
        y_pred_mean, y_pred_std = gpr.predict(X_test, return_std=True)
        df.loc[df["review"].isnull(), "pred_mean"] = y_pred_mean
        df.loc[df["review"].isnull(), "pred_std"] = y_pred_std

        members = []
        for i, row in df.iterrows():

            if isna(row["review"]):

                members.append(
                    PoseSetMember(
                        pose_id=row["pose_id"],
                        parent=self,
                        predicted_score=row["pred_mean"],
                        prediction_confidence=row["pred_std"],
                    )
                )

            else:

                members.append(
                    PoseSetMember(
                        pose_id=row["pose_id"],
                        parent=self,
                        predicted_score=None,
                        prediction_confidence=None,
                    )
                )

        PoseSetMember.objects.bulk_create(
            members,
            update_conflicts=True,
            unique_fields=["parent", "pose"],
            update_fields=["predicted_score", "prediction_confidence"],
        )

    def clear_reviews(self):
        self.poses.update(review=None)


class PoseSetMember(AbstractModel):

    class Meta:
        unique_together = ("parent", "pose")

    parent = models.ForeignKey(
        "PoseSet", on_delete=models.CASCADE, related_name="poses"
    )
    pose = models.ForeignKey("Pose", on_delete=models.PROTECT, related_name="+")

    pc1 = models.FloatField(null=True)
    pc2 = models.FloatField(null=True)

    predicted_score = models.FloatField(
        null=True,
        # validators=[
        #     MinValueValidator(0.0),
        #     MaxValueValidator(1.0)
        # ],
    )

    prediction_confidence = models.FloatField(
        null=True,
        validators=[MinValueValidator(0.0), MaxValueValidator(1.0)],
    )

    review_types = {
        0: "Bad",
        1: "Good",
    }

    review = models.IntegerField(null=True, choices=review_types)

    _style = "pose"
    _exclude_from_index = True


### ANNOTATION


class Tag(AbstractModel):
    class Meta:
        unique_together = ("name", "type")

    name = models.CharField(max_length=60, blank=False)
    type = models.ForeignKey("TagType", on_delete=models.RESTRICT, related_name="tags")

    _name_field = "name"
    _list_view_fields = ["name", "type"]
    _style = "annotation"

    def __str__(self):
        return f'"{self.name}" [{self.type.name}]'

    _field_render_types = AbstractModel._field_render_types.copy()
    _field_render_types.update(
        {
            "poses": dict(
                type=FieldRenderType.TOGGLE_CARD,
                content=ContentRenderType.INSTANCE_PILL,
                follow_related="pose",
            ),
        }
    )


class TagType(AbstractModel):

    class Meta:
        unique_together = ("name", "origin")

    name = models.CharField(max_length=30, blank=False)
    origin = models.CharField(max_length=30, blank=True, null=True)

    _style = "annotation"
    _name_field = "name"

    def __str__(self):
        return f"TT{self.id} {self.name} ({self.origin})"


class PoseTag(AbstractModel):

    class Meta:
        unique_together = ("pose", "tag")

    pose = models.ForeignKey("Pose", related_name="tags", on_delete=models.CASCADE)
    tag = models.ForeignKey("Tag", related_name="poses", on_delete=models.CASCADE)

    _exclude_from_index = True


class CompoundTag(AbstractModel):

    class Meta:
        unique_together = ("compound", "tag")

    compound = models.ForeignKey(
        "Compound", related_name="tags", on_delete=models.CASCADE
    )
    tag = models.ForeignKey("Tag", related_name="compounds", on_delete=models.CASCADE)

    _exclude_from_index = True

    # def __str__(self):
    #     return f"CT{self.id} \"{self.}\" [{self.type}]"


class Inspiration(AbstractModel):
    class Meta:
        unique_together = ("original", "derivative")

    original = models.ForeignKey("Pose", on_delete=models.CASCADE, related_name="+")
    derivative = models.ForeignKey("Pose", on_delete=models.CASCADE, related_name="+")

    _style = "annotation"


class Placement(AbstractModel):
    class Meta:
        unique_together = ("structure", "pose")

    structure = models.ForeignKey(
        "Structure", on_delete=models.RESTRICT, related_name="placements"
    )

    pose = models.OneToOneField(
        "Pose", on_delete=models.CASCADE, related_name="placement"
    )

    compound = models.ForeignKey(
        "Compound", on_delete=models.RESTRICT, related_name="placements"
    )

    # method = models.TextField(blank=True, null=True)

    metadata = models.JSONField(default=dict, blank=True)

    _name_field = "name"
    _style = "annotation"

    _list_view_fields = ["structure", "compound"]

    @property
    def name(self):
        return f"{self.structure.alias} w/ {self.compound}"


### RESOURCE MANAGEMENT


class File(AbstractModel):

    FILE_FORMATS = {
        ".pdb": "Protein Data Bank (.pdb)",
        ".mol": "Molecule (.mol)",
        ".sdf": "Structure Data File (.sdf)",
        ".csv": "Comma Separated Values (.csv)",
        ".xlsx": "Microsoft Excel (.xlsx)",
        ".json": "JavaScript Object Notation (.json)",
        ".yaml": "Yet Another Markup Language (.yaml)",
        ".html": "Hyper-Text Markup Language (.html)",
        ".smi": "SMILES (.smi)",
        ".cif": "Crystallographic Information File (.cif)",
        # zip
        # tgz
        # tar.gz
        # pkl.gz
    }

    FILE_CONTENTS = {
        "PROTEIN": "Protein structure",
        "COMPLEX": "Protein-Ligand structure",
        "LIGAND": "Ligand file",
    }

    FILE_PURPOSES = {
        "INPUT": "Input",
        "OUTPUT": "Output",
        "SCRATCH": "Scratch",
        "ARCHIVE": "Archive",
        "CONTEXT": "Context",
    }

    path = models.CharField(max_length=300, blank=False, unique=True)

    format_type = models.CharField(
        max_length=10, choices=[(k, v) for k, v in FILE_FORMATS.items()]
    )

    content_type = models.CharField(
        max_length=10, choices=[(k, v) for k, v in FILE_CONTENTS.items()]
    )
    purpose = models.CharField(
        max_length=10, choices=[(k, v) for k, v in FILE_PURPOSES.items()]
    )

    _name_field = "name"
    _style = "resources"

    _field_render_types = AbstractModel._field_render_types.copy()
    _field_render_types.update(
        {
            "content_type": dict(
                type=FieldRenderType.TABLE,
                content=ContentRenderType.TEXT_MONOSPACE,
            ),
            "path": dict(
                type=FieldRenderType.TABLE,
                content=ContentRenderType.TEXT_MONOSPACE,
                copyable=True,
            ),
            "purpose": dict(
                type=FieldRenderType.TABLE,
                content=ContentRenderType.TEXT_MONOSPACE,
            ),
            "format_type": dict(
                type=FieldRenderType.TABLE,
                content=ContentRenderType.TEXT_MONOSPACE,
            ),
        }
    )

    _list_view_fields = ["content_type", "purpose"]

    @property
    def name(self):
        return Path(self.path).name

    @property
    def as_path(self):
        return Path(self.path)


class FragalysisDownload(AbstractModel):

    from fragalysis.requests import STACKS

    STATUSES = {
        0: "PENDING",
        1: "DOWNLOADING",
        2: "LOADING",
        3: "COMPLETE",
        4: "FAILED",
    }

    target_name = models.CharField(max_length=200)
    target_access_string = models.CharField(max_length=200)

    access_token = models.CharField(max_length=32, null=True, blank=True)
    stack = models.CharField(
        max_length=60, choices={v: f"{k} ({v})" for k, v in STACKS.items()}
    )

    time_start = models.DateTimeField(auto_now_add=True)
    time_finished = models.DateTimeField(null=True)

    status = models.PositiveIntegerField(default=0, choices=STATUSES)
    message = models.TextField(null=True, blank=True)

    target = models.ForeignKey(
        "Target",
        on_delete=models.CASCADE,
        related_name="+",
        null=True,
        blank=True,
    )

    _custom_detail_view = True
    _custom_list_view = True
    _exclude_from_index = True
    _has_list_view = False

    _field_render_types = AbstractModel._field_render_types.copy()
    _field_render_types.update(
        {
            "target_name": dict(
                type=FieldRenderType.TABLE,
                content=ContentRenderType.TEXT_MONOSPACE,
                copyable=True,
            ),
            "target_access_string": dict(
                type=FieldRenderType.TABLE,
                content=ContentRenderType.TEXT_MONOSPACE,
                copyable=True,
            ),
            "stack": dict(
                type=FieldRenderType.TABLE,
                content=ContentRenderType.LINK,
            ),
            "status": dict(
                type=FieldRenderType.TABLE,
                content=ContentRenderType.TEXT_DISPLAY,
            ),
            "message": dict(
                type=FieldRenderType.CARD,
                content=ContentRenderType.TEXT_MONOSPACE,
                split="\n",
            ),
        }
    )

    # _shorthand = "FragalysisDownload"


class SdfUpload(AbstractModel):

    STATUSES = {
        0: "PENDING",
        1: "LOADING",
        2: "COMPLETE",
        3: "FAILED",
    }

    target = models.ForeignKey("Target", on_delete=models.CASCADE, related_name="+")

    protein_field_name = models.CharField(default="ref_pdb", max_length=32)
    inspirations_field_name = models.CharField(default="ref_mols", max_length=32)
    input_file = models.FileField(upload_to="media/sdf_uploads/")
    sdf_file = models.ForeignKey(
        "File", on_delete=models.RESTRICT, related_name="+", null=True
    )

    pose_set = models.OneToOneField(
        "PoseSet", on_delete=models.SET_NULL, related_name="+", null=True
    )

    compute_umap = models.BooleanField(default=True)

    time_start = models.DateTimeField(auto_now_add=True)
    time_finished = models.DateTimeField(null=True)

    status = models.PositiveIntegerField(default=0, choices=STATUSES)
    message = models.TextField(null=True, blank=True)

    pose_origins = models.CharField(
        max_length=10,
        default="COMPUTED",
        choices=[(k, v) for k, v in Pose.POSE_ORIGINS.items()],
    )

    _custom_detail_view = True
    _custom_list_view = True
    _exclude_from_index = True
    _has_list_view = False

    _field_render_types = AbstractModel._field_render_types.copy()
    _field_render_types.update(
        {
            # "target_name": dict(
            #     type=FieldRenderType.TABLE,
            #     content=ContentRenderType.TEXT_MONOSPACE,
            #     copyable=True,
            # ),
            # "target_access_string": dict(
            #     type=FieldRenderType.TABLE,
            #     content=ContentRenderType.TEXT_MONOSPACE,
            #     copyable=True,
            # ),
            "input_file": dict(
                type=FieldRenderType.TABLE,
                content=ContentRenderType.TEXT_MONOSPACE,
            ),
            "protein_field_name": dict(
                type=FieldRenderType.TABLE,
                content=ContentRenderType.TEXT_MONOSPACE,
            ),
            "inspirations_field_name": dict(
                type=FieldRenderType.TABLE,
                content=ContentRenderType.TEXT_MONOSPACE,
            ),
            "status": dict(
                type=FieldRenderType.TABLE,
                content=ContentRenderType.TEXT_DISPLAY,
            ),
            "message": dict(
                type=FieldRenderType.CARD,
                content=ContentRenderType.TEXT_MONOSPACE,
                split="\n",
            ),
            "pose_origins": dict(
                type=FieldRenderType.TABLE,
                content=ContentRenderType.TEXT_MONOSPACE,
                copyable=False,
            ),
        }
    )


MODELS = [
    Target,
    Compound,
    Pose,
    Tag,
    TagType,
    Structure,
    Placement,
    File,
    Inspiration,
    PoseSet,
    PoseSetMember,
    FragalysisDownload,
    SdfUpload,
    PoseTag,
    CompoundTag,
]
