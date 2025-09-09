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
from django.contrib.auth.models import User
from django.db.models import Avg

import sys
import mrich
import mcol
from mrich import print

from decimal import Decimal

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

    _prop_field_rendertypes = {}

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

    @property
    def poses(self):
        return Pose.objects.filter(placement__structure__target=self)


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
        from .tools import inchikey_from_smiles, sanitise_smiles, SanitisationError
        from .orm.formatters import dict_formatter

        if metadata is None:
            metadata = [None] * len(smiles)

        # get or create instance
        objects = []
        smiles_map = {}
        for s, m in zip(smiles, metadata):

            try:
                s_new = sanitise_smiles(s, verbosity=False, radical=radical)
            except SanitisationError:
                mrich.error("Could not sanitise", s)
                continue

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

    inchikey = models.CharField(max_length=27, blank=True)
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

    inspiration_distance = models.FloatField(
        validators=[MinValueValidator(0.0)], null=True
    )
    binding_energy = models.FloatField(null=True)
    ligand_energy = models.FloatField(null=False)

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
            "inspiration_distance": dict(
                type=FieldRenderType.TABLE,
                content=ContentRenderType.TEXT_MONOSPACE,
                copyable=False,
            ),
            "binding_energy": dict(
                type=FieldRenderType.TABLE,
                content=ContentRenderType.TEXT_MONOSPACE,
                copyable=False,
            ),
            "ligand_energy": dict(
                type=FieldRenderType.TABLE,
                content=ContentRenderType.TEXT_MONOSPACE,
                copyable=False,
            ),
            "tags": dict(type=FieldRenderType.HIDDEN),
            "reviews": dict(type=FieldRenderType.HIDDEN),
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

    @classmethod
    def calculate_ligand_energy(cls, mol):
        from rdkit.Chem import Mol
        from rdkit.Chem.AllChem import UFFGetMoleculeForceField, UFFOptimizeMolecule

        mol = Mol(mol)
        energy_current = UFFGetMoleculeForceField(mol).CalcEnergy()
        UFFOptimizeMolecule(mol)
        energy_minimised = UFFGetMoleculeForceField(mol).CalcEnergy()
        return energy_current - energy_minimised

    def save(self, *args, **kwargs):
        if self.ligand_energy is None and self.mol:
            self.ligand_energy = self.calculate_ligand_energy()
        super().save(*args, **kwargs)

    def __str__(self):
        tag = self.tags.filter(
            tag__type__name="Observation Code", tag__type__origin="Fragalysis"
        ).first()
        if tag:
            return f'P{self.id} "{tag.tag.name}"'
        return f"P{self.id}"


class PoseSet(AbstractModel):

    name = models.CharField(max_length=300, unique=True)

    _name_field = "name"
    _custom_detail_view = True
    _style = "pose"

    _field_render_types = AbstractModel._field_render_types.copy()
    _field_render_types.update(
        {
            "poses": dict(type=FieldRenderType.HIDDEN),
            "upload": dict(type=FieldRenderType.HIDDEN),
        }
    )

    def compute_embedding(self):

        mrich.bold(self, ".compute_embedding()")

        n = self.poses.count()

        mrich.var("#poses", n)

        if n < 3:
            mrich.warning("Not enough poses")
            return False

        with mrich.loading("imports..."):
            from math import comb
            from itertools import combinations
            from mucos import MuCOS_score
            import numpy as np
            import umap

        # get related poses
        members = self.poses.prefetch_related("pose").all()
        pose_ids = set(members.values_list("pose_id", flat=True))
        poses = {
            p.id: p
            for p in Pose.objects.filter(id__in=pose_ids).prefetch_related(
                "inspirations"
            )
        }

        n_combos = comb(n, 2)
        mrich.var("#pairs", n_combos)

        # get existing relevant pairs
        existing = set(
            PosePair.objects.filter(
                models.Q(pose1_id__in=pose_ids) | models.Q(pose2_id__in=pose_ids)
            ).values_list("pose1_id", "pose2_id")
        )
        mrich.var("#existing pairs", len(existing))

        # calculate combinations
        combos = combinations(members.values_list("pose_id", flat=True), 2)
        combos = set(tuple(sorted(t)) for t in combos) - existing

        mrich.var("#new pairs", len(combos))

        # create missing pairs
        if combos:
            new_pairs = []
            for pose1_id, pose2_id in mrich.track(combos):

                pose1 = poses[pose1_id]
                pose2 = poses[pose2_id]

                molecular_similarity = MuCOS_score(pose1.mol, pose2.mol)

                # inspiration similarity
                ref_ids1 = set(pose1.inspirations.values_list("id", flat=True))
                ref_ids2 = set(pose2.inspirations.values_list("id", flat=True))
                all_refs = ref_ids1.union(ref_ids2)
                shared_refs = ref_ids1.intersection(ref_ids2)
                if not all_refs:
                    inspiration_similarity = 0
                else:
                    inspiration_similarity = len(shared_refs) / len(all_refs)

                new_pairs.append(
                    PosePair(
                        pose1_id=pose1_id,
                        pose2_id=pose2_id,
                        molecular_similarity=molecular_similarity,
                        inspiration_similarity=inspiration_similarity,
                    )
                )

            PosePair.objects.bulk_create(new_pairs)

        # get all relevant pairs
        pairs = PosePair.objects.filter(
            models.Q(pose1_id__in=pose_ids) | models.Q(pose2_id__in=pose_ids)
        )

        # lookup dict from pose ID
        member_indices = {m.pose_id: i for i, m in enumerate(members)}

        # compute distance matrix
        distance_matrix = np.zeros((n, n))
        for pair in pairs:
            i = member_indices.get(pair.pose1_id)
            j = member_indices.get(pair.pose2_id)

            if not i:
                mrich.error(f"Could not get index for Pose w/ id={pair.pose1_id}")
            if not j:
                mrich.error(f"Could not get index for Pose w/ id={pair.pose2_id}")

            dist = pair.distance
            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist

        """ 
        ## ADD TO SIMILARITY COMPARISON
            - distance_score similarity
            - energy_score similarity
        """

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

    @property
    def unreviewed(self):
        return self.poses.filter(pose__reviews__isnull=True)

    def generate_umap_fig(
        self, force_embedding: bool = False, force_extrapolation: bool = False
    ):

        import plotly.graph_objects as go
        from pandas import DataFrame
        import numpy as np

        # Compute UMAP embedding if needed
        if (
            self.poses.filter(pc1__isnull=False).count() != self.poses.count()
            or force_embedding
        ):
            self.compute_embedding()
        else:
            mrich.debug("No UMAP embedding needed")

        members = self.poses

        annotated = members.annotate(review=Avg("pose__reviews__review"))

        mrich.var(
            "reviewed w/ prediction",
            annotated.filter(
                review__isnull=False, predicted_score__isnull=False
            ).count(),
        )
        mrich.var(
            "unreviewed w/o prediction",
            annotated.filter(review__isnull=True, predicted_score__isnull=True).count(),
        )

        if (
            annotated.filter(
                review__isnull=False, predicted_score__isnull=False
            ).count()
            or annotated.filter(
                review__isnull=True, predicted_score__isnull=True
            ).count()
            or force_extrapolation
        ):
            self.extrapolate_reviews()
        else:
            mrich.debug("No extrapolation needed")

        data = annotated.values(
            "pose_id",
            "pc1",
            "pc2",
            "predicted_score",
            "prediction_uncertainty",
            "review",
        )  # .annotate(review=Avg("pose__reviews__review"))

        df = DataFrame(data)

        fig = go.Figure()

        if members.count():

            df_train = df[~df["review"].isnull()]
            df_test = df[df["review"].isnull()]

            if not len(df_train):
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

                # HUMAN REVIEWS

                customdata = [
                    (pose_id, review, 1)
                    for pose_id, review in df_train[["pose_id", "review"]].values
                ]

                fig.add_trace(
                    go.Scatter(
                        name="Reviewed",
                        x=df_train["pc1"],
                        y=df_train["pc2"],
                        mode="markers",
                        marker=dict(
                            color=(1 + df_train["review"]) / 2,
                            symbol="star",
                            size=10,
                            colorscale=[
                                [0.0, "red"],
                                [0.5, "white"],
                                [1.0, "rgb(0,255,0)"],
                            ],
                            line=dict(color="black", width=1),
                            cmin=0,
                            cmax=1,
                        ),
                        hovertemplate="P%{customdata[0]} [%{customdata[1]}]<br>PC1: %{x}<br>PC2: %{y}<extra></extra>",
                        customdata=customdata,
                    )
                )

                # EXTRAPOLATED

                data = df_test[
                    ["pose_id", "predicted_score", "prediction_uncertainty"]
                ].values

                customdata = [(int(a), b, c) for a, b, c in data]

                fig.add_trace(
                    go.Scatter(
                        name="Predicted",
                        x=df_test["pc1"],
                        y=df_test["pc2"],
                        mode="markers",
                        marker=dict(
                            color=(1 + df_test["predicted_score"]) / 2,
                            colorscale=[
                                [0.0, "red"],
                                [0.5, "white"],
                                [1.0, "rgb(0,255,0)"],
                            ],
                            opacity=np.maximum(
                                1 - df_test["prediction_uncertainty"], 0.2
                            ),
                            line=dict(color="black", width=1),
                            cmin=0,
                            cmax=1,
                        ),
                        hovertemplate="P%{customdata[0]} [%{customdata[1]:.2f} +/- %{customdata[2]:.2f}]<br>PC1: %{x}<br>PC2: %{y}<br>Predicted: %{customdata[1]}<br>Uncertainty: %{customdata[2]}<extra></extra>",
                        customdata=customdata,
                    )
                )

        fig.update_layout(
            margin=dict(l=20, r=20, t=20, b=20),
            paper_bgcolor="rgba(0,0,0,0)",
        )

        return fig

    # def generate_score_

    def extrapolate_reviews(self, debug: bool = False):

        if not self.poses.count():
            return

        mrich.bold(self, ".extrapolate_reviews()")

        from pandas import DataFrame, isna

        data = (
            self.poses.select_related("pose")
            .values(
                "pose_id",
                "pc1",
                "pc2",
                "pose__inspiration_distance",
                "pose__binding_energy",
                "pose__ligand_energy",
            )
            .annotate(review=Avg("pose__reviews__review"))
        )

        df = DataFrame(data)

        df["pose__inspiration_distance"] = df["pose__inspiration_distance"].fillna(10)
        df["pose__binding_energy"] = df["pose__inspiration_distance"].fillna(100)
        df["pose__ligand_energy"] = df["pose__ligand_energy"].fillna(200)

        df_train = df[df["review"].notnull()]
        df_test = df[df["review"].isnull()]

        if debug:
            mrich.var("#total", len(df))
            mrich.var("#train", len(df_train))
            mrich.var("#test", len(df_test))

        if not len(df_train):
            mrich.warning("Can't extrapolate, no reviews")
            self.clear_predictions()
            return

        if not len(df_test):
            mrich.warning("Can't extrapolate, no unreviewed")
            self.clear_predictions()
            return

        cols = [
            "pc1",
            "pc2",
            "pose__inspiration_distance",
            "pose__binding_energy",
            "pose__ligand_energy",
        ]

        X_train = df_train[cols].values
        y_train = df_train["review"].values
        X_test = df_test[cols].values

        with mrich.loading("imports..."):
            from sklearn.gaussian_process import GaussianProcessRegressor
            from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C
            import numpy as np

        kernel = C(1.0, (1e-3, 1e3)) * RBF(
            length_scale=1.0, length_scale_bounds=(1e-2, 1e2)
        )
        gpr = GaussianProcessRegressor(
            kernel=kernel, n_restarts_optimizer=10, alpha=1e-2
        )
        gpr.fit(X_train, y_train)

        y_pred_mean, y_pred_std = gpr.predict(X_test, return_std=True)

        if debug:
            mrich.var("mean prediction STD", np.mean(y_pred_std))
            mrich.var("max prediction STD", np.max(y_pred_std))
            mrich.var("min prediction STD", np.min(y_pred_std))

        df.loc[df["review"].isnull(), "pred_mean"] = [
            max(-1, min(1, v)) for v in y_pred_mean
        ]
        df.loc[df["review"].isnull(), "pred_std"] = y_pred_std

        members = []
        for i, row in df.iterrows():

            if isna(row["review"]):

                members.append(
                    PoseSetMember(
                        pose_id=row["pose_id"],
                        parent=self,
                        predicted_score=row["pred_mean"],
                        prediction_uncertainty=row["pred_std"],
                    )
                )

            else:

                members.append(
                    PoseSetMember(
                        pose_id=row["pose_id"],
                        parent=self,
                        predicted_score=None,
                        prediction_uncertainty=None,
                    )
                )

        PoseSetMember.objects.bulk_create(
            members,
            update_conflicts=True,
            unique_fields=["parent", "pose"],
            update_fields=["predicted_score", "prediction_uncertainty"],
        )

    def clear_predictions(self):
        self.poses.update(predicted_score=None, prediction_uncertainty=None)

    def clear_embedding(self):
        self.poses.update(pc1=None, pc2=None)


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
        validators=[MinValueValidator(-1.0), MaxValueValidator(1.0)],
    )

    prediction_uncertainty = models.FloatField(
        null=True,
        validators=[MinValueValidator(0.0), MaxValueValidator(1.0)],
    )

    _style = "pose"
    _exclude_from_index = True


class PoseReview(AbstractModel):

    class Meta:
        unique_together = ("user", "pose")

    user = models.ForeignKey(User, on_delete=models.CASCADE, related_name="+")
    pose = models.ForeignKey(Pose, on_delete=models.CASCADE, related_name="reviews")

    REVIEW_TYPES = {
        Decimal("-1.00"): "Very Bad",
        Decimal("-0.75"): "Bad",
        Decimal("0.00"): "Neutral",
        Decimal("0.25"): "Interesting",
        Decimal("0.75"): "Good",
        Decimal("1.00"): "Great",
    }

    review = models.DecimalField(
        max_digits=4,
        decimal_places=2,
        choices=[(k, v) for k, v in REVIEW_TYPES.items()],
    )

    _style = "annotation"
    _exclude_from_index = True

    def __str__(self):
        return (
            f"PR{self.id} {self.user}: P{self.pose.id} => {self.get_review_display()}"
        )


class PosePair(AbstractModel):

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=["pose1", "pose2"], name="unique_pair"),
            models.CheckConstraint(
                check=models.Q(pose1_id__lt=models.F("pose2_id")), name="pose_order"
            ),
        ]

    pose1 = models.ForeignKey("Pose", related_name="+", on_delete=models.CASCADE)
    pose2 = models.ForeignKey("Pose", related_name="+", on_delete=models.CASCADE)

    molecular_similarity = models.FloatField(
        validators=[MinValueValidator(0.0), MaxValueValidator(1.0)]
    )
    inspiration_similarity = models.FloatField(
        validators=[MinValueValidator(0.0), MaxValueValidator(1.0)]
    )

    _exclude_from_index = True

    @property
    def similarity(self):
        return 0.8 * self.molecular_similarity + 0.2 * self.inspiration_similarity

    @property
    def distance(self):
        return 1 - self.similarity

    def __str__(self):
        return f"PP{self.id} (P{self.pose1_id}, P{self.pose2_id})"


### ANNOTATION


class Tag(AbstractModel):
    class Meta:
        unique_together = ("name", "type")

    name = models.CharField(max_length=60, blank=False)
    type = models.ForeignKey("TagType", on_delete=models.RESTRICT, related_name="tags")

    _name_field = "name"
    _list_view_fields = ["name", "type"]
    _style = "annotation"
    _exclude_from_index = True

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
    _exclude_from_index = True

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
    _exclude_from_index = True


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


### PROCUREMENT


class Supplier(AbstractModel):

    name = models.CharField(max_length=30, unique=True)

    _shorthand = "Su"
    _style = "quoting"
    _name_field = "name"


class Quote(AbstractModel):

    class Meta:
        unique_together = ("amount", "supplier", "catalogue", "entry")

    supplier = models.ForeignKey(
        "Supplier", on_delete=models.CASCADE, related_name="quotes"
    )
    compound = models.ForeignKey(
        "Compound", on_delete=models.CASCADE, related_name="quotes"
    )

    amount = models.FloatField()  # mg
    catalogue = models.CharField(max_length=90, null=True, blank=True)
    entry = models.CharField(max_length=60)
    lead_time = models.FloatField(null=True, blank=True)
    date = models.DateField()
    purity = models.FloatField(null=True, blank=True)

    price = models.DecimalField(
        max_digits=8,
        decimal_places=2,
    )

    currency = models.CharField(max_length=3, choices=CURRENCIES)

    smiles = models.CharField(max_length=300, unique=False)

    mol = models.GeneratedField(
        expression=MolFromSmiles("smiles", "mol"),
        output_field=MolField(blank=True),
        db_persist=True,
    )

    _style = "quoting"

    _list_view_fields = ["supplier", "compound", "entry"]

    _field_render_types = AbstractModel._field_render_types.copy()
    _field_render_types.update(
        {
            "compound": dict(
                type=FieldRenderType.TABLE,
                content=ContentRenderType.INSTANCE_PILL,
                order=-1,
            ),
            "supplier": dict(
                type=FieldRenderType.TABLE,
                content=ContentRenderType.INSTANCE_PILL,
                order=-1,
            ),
            "entry": dict(
                type=FieldRenderType.TABLE,
                content=ContentRenderType.TEXT_MONOSPACE,
                copyable=True,
            ),
            "lead_time": dict(
                type=FieldRenderType.TABLE,
                content=ContentRenderType.TEXT_MONOSPACE,
                copyable=False,
            ),
            "currency": dict(
                type=FieldRenderType.HIDDEN,
            ),
            "price": dict(
                type=FieldRenderType.HIDDEN,
            ),
            "mol": dict(
                type=FieldRenderType.HIDDEN,
            ),
        }
    )

    _prop_field_rendertypes = {
        "price": dict(
            prop="price_str",
            type=FieldRenderType.TABLE,
            content=ContentRenderType.TEXT_MONOSPACE,
            copyable=True,
        ),
    }

    @property
    def price_str(self) -> str:
        return f"{self.get_currency_display()}{self.price}"

    def __str__(self):
        """Unformatted string representation"""
        if self.purity:
            purity = f" @ {self.purity:.0%}"
        else:
            purity = ""

        if self.supplier == "Stock":
            s = f"Q{self.id} {self.compound} In Stock: {self.amount:}mg{purity}"
        else:
            s = f"Q{self.id} {self.compound} {self.entry} {self.amount}mg{purity} = {self.currency}{self.price}"

        if self.lead_time:
            s += f" ({self.lead_time} days)"

        return s


### SCORES


class CompoundScore(AbstractModel):

    compound = models.ForeignKey(
        "Compound", on_delete=models.CASCADE, related_name="scores"
    )
    type = models.ForeignKey(
        "ScoreType", on_delete=models.CASCADE, related_name="compound_scores"
    )
    target = models.ForeignKey(
        "Target", on_delete=models.CASCADE, related_name="compound_scores"
    )
    value = models.FloatField()
    unit = models.CharField(max_length=10)

    _style = "annotation"

    _field_render_types = AbstractModel._field_render_types.copy()
    _field_render_types.update(
        {
            "unit": dict(
                type=FieldRenderType.TABLE,
                content=ContentRenderType.TEXT_MONOSPACE,
                copyable=False,
            ),
        }
    )

    _list_view_fields = ["compound", "type", "target", "value", "unit"]

    def __str__(self):
        return f"CS{self.id} ST{self.type_id} {self.value} {self.unit} [C{self.compound_id}]"


class PoseScore(AbstractModel):

    pose = models.ForeignKey("Pose", on_delete=models.CASCADE, related_name="scores")
    type = models.ForeignKey(
        "ScoreType", on_delete=models.CASCADE, related_name="pose_scores"
    )
    value = models.FloatField()
    unit = models.CharField(max_length=10)

    _style = "annotation"

    _field_render_types = AbstractModel._field_render_types.copy()
    _field_render_types.update(
        {
            "unit": dict(
                type=FieldRenderType.TABLE,
                content=ContentRenderType.TEXT_MONOSPACE,
                copyable=False,
            ),
        }
    )

    _list_view_fields = ["pose", "type", "value", "unit"]

    def __str__(self):
        return (
            f"PS{self.id} ST{self.type_id} {self.value} {self.unit} [P{self.pose_id}]"
        )


class ScoreType(AbstractModel):

    SCORE_ORIGINS = {
        "EXPERIMENT": "Experimental measurement",
        "COMPUTED": "Virtual / computed measurement",
        "MANUAL": "Human-assigned value",
    }

    name = models.CharField(max_length=60, unique=True)
    type = models.CharField(max_length=10, choices=SCORE_ORIGINS)

    _style = "annotation"

    def __str__(self):
        return f"ST{self.id} '{self.name}' [{self.type}]"


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
        "META": "Metadata",
        "PLOT": "Graph",
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

    from fragalysis.requests.urls import STACKS

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
    binding_energy_field_name = models.CharField(
        default="energy_score", max_length=32, blank=True, null=True
    )
    inspiration_distance_field_name = models.CharField(
        default="distance_score", max_length=32, blank=True, null=True
    )
    # ligand_energy_field_name = models.CharField(default="ref_mols", max_length=32, blank=True, null=True)
    input_file = models.FileField(upload_to="media/sdf_uploads/")
    sdf_file = models.ForeignKey(
        "File", on_delete=models.RESTRICT, related_name="+", null=True
    )

    pose_set = models.OneToOneField(
        "PoseSet", on_delete=models.SET_NULL, related_name="upload", null=True
    )

    compute_embedding = models.BooleanField(default=True)

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


class ScoreUpload(AbstractModel):

    STATUSES = {
        0: "PENDING",
        1: "LOADING",
        2: "COMPLETE",
        3: "FAILED",
    }

    SCORE_TYPE = {
        0: "Compound",
        1: "Pose",
    }

    target = models.ForeignKey("Target", on_delete=models.CASCADE, related_name="+")
    score_type = models.ForeignKey(
        "ScoreType", on_delete=models.CASCADE, related_name="+"
    )

    related_model = models.PositiveIntegerField(default=0, choices=SCORE_TYPE)
    identifier_column_name = models.CharField(max_length=32)
    separator = models.CharField(max_length=3, default=",")
    value_field = models.CharField(max_length=120, blank=False, null=False)
    unit_field = models.CharField(max_length=120, blank=False, null=False)

    input_file = models.FileField(upload_to="media/score_uploads/")

    file = models.ForeignKey(
        "File", on_delete=models.RESTRICT, related_name="+", null=True
    )

    time_start = models.DateTimeField(auto_now_add=True)
    time_finished = models.DateTimeField(null=True)

    status = models.PositiveIntegerField(default=0, choices=STATUSES)
    message = models.TextField(null=True, blank=True)

    _custom_detail_view = True
    _custom_list_view = True
    _exclude_from_index = True
    _has_list_view = False

    # _field_render_types = AbstractModel._field_render_types.copy()
    # _field_render_types.update(
    #     {
    #         # "target_name": dict(
    #         #     type=FieldRenderType.TABLE,
    #         #     content=ContentRenderType.TEXT_MONOSPACE,
    #         #     copyable=True,
    #         # ),
    #         # "target_access_string": dict(
    #         #     type=FieldRenderType.TABLE,
    #         #     content=ContentRenderType.TEXT_MONOSPACE,
    #         #     copyable=True,
    #         # ),
    #         "input_file": dict(
    #             type=FieldRenderType.TABLE,
    #             content=ContentRenderType.TEXT_MONOSPACE,
    #         ),
    #         "protein_field_name": dict(
    #             type=FieldRenderType.TABLE,
    #             content=ContentRenderType.TEXT_MONOSPACE,
    #         ),
    #         "inspirations_field_name": dict(
    #             type=FieldRenderType.TABLE,
    #             content=ContentRenderType.TEXT_MONOSPACE,
    #         ),
    #         "status": dict(
    #             type=FieldRenderType.TABLE,
    #             content=ContentRenderType.TEXT_DISPLAY,
    #         ),
    #         "message": dict(
    #             type=FieldRenderType.CARD,
    #             content=ContentRenderType.TEXT_MONOSPACE,
    #             split="\n",
    #         ),
    #         "pose_origins": dict(
    #             type=FieldRenderType.TABLE,
    #             content=ContentRenderType.TEXT_MONOSPACE,
    #             copyable=False,
    #         ),
    #     }
    # )


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
    ScoreUpload,
    PoseTag,
    CompoundTag,
    PoseReview,
    PosePair,
    Supplier,
    Quote,
    CompoundScore,
    PoseScore,
    ScoreType,
]
