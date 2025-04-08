"""This file defines the base models for the HIPPO Database schema using django, these classes are then subclassed in other files to create Python wrappers"""

from pathlib import Path
from datetime import date

from django.db import models
from django.core.validators import MinValueValidator, MaxValueValidator

from .orm.expressions import MolFromSmiles, MolPatternBfpFromSmiles
from .orm.fields import MolField
from .orm.validators import validate_list_of_integers, validate_coord

from molparse.rdkit.features import FEATURE_FAMILIES, INTERACTION_TYPES

import sys
import mrich

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

    _custom_detail_view = False

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
        "inchikey": dict(
            type=FieldRenderType.TABLE,
            content=ContentRenderType.TEXT_MONOSPACE,
            copyable=True,
        ),
        "mol": dict(type=FieldRenderType.TABLE, content=ContentRenderType.MOL_2D_SVG),
        "pdbblock": dict(type=FieldRenderType.HIDDEN),
        "metadata": dict(
            type=FieldRenderType.TOGGLE_CARD, content=ContentRenderType.DICT_TABLE
        ),
    }

    _list_view_fields = ["smiles", "inchikey", "mol"]

    @property
    def model_pill_html(self):
        return f"""<div class="model-pill" style="
        background:var(--color-{self.__name__.lower()}, 
                   var(--color-{self._parent_module}));
        color:var(--text-color-{self.__name__.lower()}, 
              var(--text-color-{self._parent_module}));
        ">{self}</div>"""

    def get_absolute_url(self):
        return reverse(f"{self.__name__.lower()}_detail", args=[str(self.id)])

    def summary(self):

        from rich.panel import Panel
        from rich.box import SIMPLE_HEAVY
        from rich.table import Table

        fields = self.get_wrapped_field_names()

        table = Table(title=self.__rich__(), box=SIMPLE_HEAVY)
        table.add_column("Field", style="var_name")
        table.add_column("Value", style="result")

        table.add_row(f"[bold]Model", f"[bold var_type]{self.__name__}")

        for field in fields:
            value = getattr(self, field)

            if not value:
                continue

            if hasattr(value, "__rich__"):
                s = value.__rich__()
            else:
                s = str(value)
            table.add_row(f"[bold]{field}", s)

        panel = Panel(table, expand=False)
        mrich.print(panel)

    def clean_and_save(self, debug: bool = False):
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

        if sh := self.shorthand:
            s = f"{self.shorthand}{id_str}"
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

    name = models.CharField(max_length=60, blank=True, unique=True)
    metadata = models.JSONField(default=dict, blank=True)

    _shorthand = "T"
    _name_field = "name"


class Structure(AbstractModel):

    alias = models.CharField(max_length=60, blank=True, unique=True, null=True)
    pdbblock = models.TextField(blank=True, null=True)

    files = models.ManyToManyField("File", related_name="_structures")

    target = models.ForeignKey(
        "Target", on_delete=models.CASCADE, related_name="_structures"
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

    poses = models.ManyToManyField(
        "Pose",
        through="Placement",
        through_fields=("structure", "pose"),
        related_name="_structures",
    )

    metadata = models.JSONField(default=dict, blank=True)

    _shorthand = "S"
    _name_field = "alias"


### LIGAND (2D)


class Compound(AbstractModel):

    inchikey = models.CharField(max_length=27, unique=True)
    alias = models.CharField(max_length=60, blank=True, unique=True, null=True)
    smiles = models.CharField(max_length=300, unique=True)
    metadata = models.JSONField(default=dict, blank=True)

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

    scaffolds = models.ManyToManyField(
        "Compound", related_name="_elaborations", blank=True
    )

    tags = models.ManyToManyField("Tag", related_name="_compounds", blank=True)

    _shorthand = "C"
    _name_field = "alias"

    _field_render_types = AbstractModel._field_render_types.copy()
    _field_render_types.update(
        {
            "pattern_bfp": dict(type=FieldRenderType.HIDDEN),
            # "mol": dict(type=FieldRenderType.HIDDEN),
        }
    )

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


### LIGAND (3D)


class Pose(AbstractModel):

    inchikey = models.CharField(max_length=27, blank=True)
    alias = models.CharField(max_length=60, blank=True, unique=True, null=True)
    smiles = models.CharField(max_length=300, blank=True, null=True)

    files = models.ManyToManyField("File", related_name="_poses")

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
        related_name="_derivatives",
    )
    tags = models.ManyToManyField("Tag", related_name="_poses")

    _shorthand = "P"
    _name_field = "alias"
    _custom_detail_view = True

    _field_render_types = AbstractModel._field_render_types.copy()
    _field_render_types.update(
        {
            "origin": dict(
                type=FieldRenderType.TABLE,
                content=ContentRenderType.TEXT_MONOSPACE,
                copyable=False,
            ),
            "_structures": dict(
                type=FieldRenderType.TABLE,
                content=ContentRenderType.INSTANCE_PILL,
                zero_index=True,
            ),
            "_placements": dict(
                type=FieldRenderType.TABLE,
                content=ContentRenderType.INSTANCE_PILL,
                zero_index=True,
            ),
        }
    )

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

    name = models.CharField(max_length=60, unique=True)

    _name_field = "name"
    _custom_detail_view = True

    def calculate_pca(self):

        with mrich.loading("imports..."):
            from ..tools import get_cfps
            from sklearn.decomposition import PCA
            import numpy as np
            import pandas as pd

        members = list(self.poses.values("id", "pose"))

        pose_ids = [d["pose"] for d in members]

        poses = Pose.objects.filter(id__in=pose_ids)

        pose_mols = {d["id"]: d["mol"] for d in poses.values("id", "mol")}

        for member in members:
            mol = pose_mols[member["pose"]]
            member["mol"] = mol

        df = pd.DataFrame(members)

        with mrich.loading("Getting Compound fingerprints"):
            df["FP"] = df["mol"].map(get_cfps)

        X = np.array([x.fp for x in df["FP"]])

        with mrich.loading("Computing PCA"):
            pca = PCA(n_components=2, random_state=0)
            pca_fit = pca.fit_transform(X)

        df["PC1"] = pca_fit.T[0]
        df["PC2"] = pca_fit.T[1]

        members = []

        for i, row in df.iterrows():

            members.append(
                PoseSetMember(
                    id=row["id"],
                    pc1=row["PC1"],
                    pc2=row["PC2"],
                )
            )

        PoseSetMember.objects.bulk_update(members, fields=["pc1", "pc2"])

    def generate_tsnee_fig(self):

        import plotly.graph_objects as go

        n_valid = [v for v in self.poses.values_list("pc1", flat=True) if v is not None]

        if n_valid != self.poses.count():
            self.calculate_pca()

        members = self.poses

        fig = go.Figure()

        x = [m.pc1 for m in members]
        y = [m.pc2 for m in members]

        text = []
        colors = []
        pose_ids = []
        for m in members:

            pose = m.pose

            text.append(str(pose))
            pose_ids.append(pose.id)

            match m.review:
                case None:
                    colors.append("gray")
                case "GOOD":
                    colors.append("rgb(0,255,0)")
                case "BAD":
                    colors.append("red")

        trace = go.Scatter(
            x=x,
            y=y,
            mode="markers",
            marker=dict(
                size=12,
                color=colors,
            ),
            hovertemplate="%{text}<br>PC1: %{x}<br>PC2: %{y}<extra></extra>",
            text=text,
            customdata=pose_ids,
        )

        fig.add_trace(trace)

        fig.update_layout(
            margin=dict(l=20, r=20, t=20, b=20),
            paper_bgcolor="rgba(0,0,0,0)",
        )

        return fig


class PoseSetMember(AbstractModel):

    parent = models.ForeignKey(
        "PoseSet", on_delete=models.CASCADE, related_name="_poses"
    )
    pose = models.ForeignKey("Pose", on_delete=models.PROTECT, related_name="+")

    pc1 = models.FloatField(null=True)
    pc2 = models.FloatField(null=True)

    review_types = {
        "BAD": "Bad",
        "GOOD": "Good",
    }

    review = models.CharField(max_length=4, null=True, choices=review_types)


### ANNOTATION


class Tag(AbstractModel):
    class Meta:
        unique_together = ("name", "type")

    name = models.CharField(max_length=60, blank=False)
    type = models.ForeignKey("TagType", on_delete=models.RESTRICT, related_name="_tags")

    _name_field = "name"
    _list_view_fields = ["name", "type"]

    def __str__(self):
        return f'"{self.name}" [{self.type.name}]'


class TagType(AbstractModel):

    name = models.CharField(max_length=30, blank=False, unique=True)
    origin = models.CharField(max_length=30, blank=True, null=True)

    _name_field = "name"


class Inspiration(AbstractModel):
    class Meta:
        unique_together = ("original", "derivative")

    original = models.ForeignKey("Pose", on_delete=models.CASCADE, related_name="+")
    derivative = models.ForeignKey("Pose", on_delete=models.CASCADE, related_name="+")


class Placement(AbstractModel):
    class Meta:
        unique_together = ("structure", "pose")

    structure = models.ForeignKey(
        "Structure", on_delete=models.RESTRICT, related_name="_placements"
    )

    pose = models.ForeignKey(
        "Pose", on_delete=models.CASCADE, related_name="_placements"
    )

    compound = models.ForeignKey(
        "Compound", on_delete=models.RESTRICT, related_name="_placements"
    )

    method = models.TextField(blank=True, null=True)

    metadata = models.JSONField(default=dict, blank=True)


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

    @property
    def name(self):
        return Path(self.path).name

    @property
    def as_path(self):
        return Path(self.path)


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
]
