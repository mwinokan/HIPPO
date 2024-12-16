"""This file defines the base models for the HIPPO Database schema using django, these classes are then subclassed in other files to create Python wrappers"""

from pathlib import Path
from datetime import date

from django.db import models
from django.core.validators import MinValueValidator, MaxValueValidator

from .orm.expressions import MolFromSmiles, MolPatternBfpFromSmiles
from .orm.fields import MolField
from .orm.validators import validate_list_of_integers, validate_coord
from .abstract import AbstractModel

from molparse.rdkit.features import FEATURE_FAMILIES, INTERACTION_TYPES

CURRENCIES = {
    "EUR": "€",
    "USD": "$",
    "GBP": "£",
}

### PROTEIN


class TargetModel(AbstractModel):
    class Meta:
        app_label = "hippo"
        abstract = True

    name = models.CharField(max_length=60, blank=True, unique=True)
    metadata = models.JSONField(default=dict, blank=True)

    _shorthand = "T"
    _name_field = "name"


class StructureModel(AbstractModel):
    class Meta:
        app_label = "hippo"
        abstract = True

    alias = models.CharField(max_length=60, blank=True, unique=True, null=True)
    pdbblock = models.TextField(blank=True, null=True)
    # path = models.FilePathField(Path("/"), max_length=200, unique=True)

    _files = models.ManyToManyField("File", related_name="_structures")

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

    _poses = models.ManyToManyField(
        "Pose",
        through="Placement",
        through_fields=("structure", "pose"),
        related_name="_structures",
    )

    metadata = models.JSONField(default=dict, blank=True)

    _shorthand = "S"
    _name_field = "alias"


### LIGAND (2D)


class CompoundModel(AbstractModel):
    class Meta:
        app_label = "hippo"
        abstract = True

    inchikey = models.CharField(max_length=27, unique=True)
    alias = models.CharField(max_length=60, blank=True, unique=True, null=True)
    smiles = models.CharField(max_length=300)
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

    _suppliers = models.ManyToManyField(
        "Supplier",
        through="Quote",
        through_fields=("compound", "supplier"),
        related_name="_compounds",
    )

    _scaffolds = models.ManyToManyField(
        "Compound", related_name="_elaborations", blank=True
    )
    _tags = models.ManyToManyField("Tag", related_name="_compounds", blank=True)

    _score_types = models.ManyToManyField(
        "CompoundScoreType",
        through="CompoundScore",
        through_fields=("compound", "score_type"),
        related_name="_compounds",
    )

    _shorthand = "C"
    _name_field = "alias"


class CompoundScoreTypeModel(AbstractModel):
    class Meta:
        abstract = True

    name = models.CharField(max_length=30)
    method = models.CharField(max_length=30)
    unit = models.CharField(max_length=30)
    description = models.CharField(max_length=120)


class CompoundScoreModel(AbstractModel):
    class Meta:
        abstract = True

    compound = models.ForeignKey(
        "Compound", on_delete=models.CASCADE, related_name="_scores"
    )
    score_type = models.ForeignKey(
        "CompoundScoreType", on_delete=models.CASCADE, related_name="_values"
    )

    value = models.FloatField()


### LIGAND (3D)


class PoseModel(AbstractModel):
    class Meta:
        app_label = "hippo"
        abstract = True

    inchikey = models.CharField(max_length=27, blank=True)
    alias = models.CharField(max_length=60, blank=True, unique=True, null=True)
    smiles = models.CharField(max_length=300, blank=True, null=True)

    # reference = models.ForeignKey(
    #     "Pose", on_delete=models.SET_NULL, blank=True, null=True, related_name="+"
    # )

    # structure = models.ForeignKey("Structure", on_delete=models.CASCADE, blank=False, null=False, related_name="_poses")

    _files = models.ManyToManyField("File", related_name="_poses")

    # mol_path = models.FilePathField(Path("/"), max_length=200, unique=True)
    # complex_path = models.FilePathField(Path("/"), max_length=200, unique=True)

    # compound = models.ForeignKey(
    #     "Compound", on_delete=models.CASCADE, related_name="_poses"
    # )
    # target = models.ForeignKey(
    #     "Target", on_delete=models.CASCADE, related_name="_poses"
    # )
    mol = MolField(blank=True)
    # energy_score = models.FloatField(blank=True, null=True)
    # distance_score = models.FloatField(blank=True, null=True)
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

    _score_types = models.ManyToManyField(
        "PoseScoreType",
        through="PoseScore",
        through_fields=("pose", "score_type"),
        related_name="_poses",
    )

    _inspirations = models.ManyToManyField(
        "Pose",
        through="Inspiration",
        through_fields=("original", "derivative"),
        related_name="_derivatives",
    )
    _tags = models.ManyToManyField("Tag", related_name="_poses")

    _shorthand = "P"
    _name_field = "alias"


class PoseScoreTypeModel(AbstractModel):
    class Meta:
        abstract = True

    name = models.CharField(max_length=30)
    method = models.CharField(max_length=30)
    unit = models.CharField(max_length=30)
    description = models.CharField(max_length=120)


class PoseScoreModel(AbstractModel):

    class Meta:
        abstract = True

    pose = models.ForeignKey("Pose", on_delete=models.CASCADE, related_name="_scores")
    score_type = models.ForeignKey(
        "PoseScoreType", on_delete=models.CASCADE, related_name="_values"
    )

    value = models.FloatField()


### INTERACTIONS


class FeatureModel(AbstractModel):

    class Meta:
        app_label = "hippo"
        abstract = True
        unique_together = (
            "family",
            "target",
            "chain_name",
            "residue_name",
            "residue_number",
            "atom_names",
        )

    family = models.CharField(max_length=30, choices={f: f for f in FEATURE_FAMILIES})

    target = models.ForeignKey(
        "Target",
        on_delete=models.CASCADE,
        related_name="_features",
    )

    chain_name = models.CharField(max_length=5)
    residue_name = models.CharField(max_length=10)
    residue_number = models.SmallIntegerField()
    atom_names = models.JSONField(default=dict)  # TODO: Validate as list

    _shorthand = "F"
    _name_field = "family"


class InteractionModel(AbstractModel):
    class Meta:
        app_label = "hippo"
        abstract = True
        unique_together = ("feature", "pose", "type", "family")

    family = models.CharField(max_length=30, choices={f: f for f in FEATURE_FAMILIES})
    type = models.CharField(
        max_length=30, choices={f: f for f in INTERACTION_TYPES.values()}
    )

    pose = models.ForeignKey(
        "Pose", on_delete=models.CASCADE, related_name="_interactions"
    )

    structure = models.ForeignKey(
        "Structure", on_delete=models.CASCADE, related_name="_interactions"
    )

    feature = models.ForeignKey(
        "Feature", on_delete=models.RESTRICT, related_name="_interactions"
    )

    atom_ids = models.JSONField(default=dict, validators=[validate_list_of_integers])
    prot_coord = models.JSONField(default=dict, validators=[validate_coord])
    lig_coord = models.JSONField(default=dict, validators=[validate_coord])
    distance = models.FloatField(validators=[MinValueValidator(0.0)])
    angle = models.FloatField(
        blank=True,
        null=True,
        validators=[MinValueValidator(0.0), MaxValueValidator(360)],
    )
    energy = models.FloatField(blank=True, null=True)

    _shorthand = "I"
    _name_field = "type"


### PROCUREMENT


class QuoteModel(AbstractModel):
    class Meta:
        app_label = "hippo"
        abstract = True
        unique_together = ("amount", "supplier", "catalogue", "entry")

    amount = models.FloatField(validators=[MinValueValidator(0.0)])
    catalogue = models.CharField(max_length=60, blank=True)
    entry = models.CharField(max_length=60, blank=True)
    lead_time = models.FloatField(
        validators=[MinValueValidator(0.0)], blank=True, null=True
    )
    price = models.DecimalField(
        validators=[MinValueValidator(0.0)],
        blank=True,
        null=True,
        max_digits=8,
        decimal_places=2,
    )
    currency = models.CharField(max_length=3, blank=True, choices=CURRENCIES)
    date = models.DateField(default=date.today)

    compound = models.ForeignKey(
        "Compound", on_delete=models.RESTRICT, related_name="_quotes"
    )

    supplier = models.ForeignKey(
        "Supplier", on_delete=models.RESTRICT, related_name="_quotes"
    )

    purity = models.FloatField(
        blank=True,
        null=True,
        validators=[MinValueValidator(0.0), MaxValueValidator(1.0)],
    )

    _shorthand = "Q"


class SupplierModel(AbstractModel):
    class Meta:
        app_label = "hippo"
        abstract = True

    name = models.CharField(max_length=30, unique=True)

    _name_field = "name"


### CHEMISTRY


class ReactionModel(AbstractModel):
    class Meta:
        app_label = "hippo"
        abstract = True

    type = models.CharField(max_length=60, blank=False)

    _product_compounds = models.ManyToManyField(
        "Compound",
        through="Product",
        through_fields=("reaction", "compound"),
        related_name="_product_reactions",
    )

    yield_fraction = models.FloatField(
        default=1.0, validators=[MinValueValidator(0.0), MaxValueValidator(1.0)]
    )

    _reactant_compounds = models.ManyToManyField(
        "Compound",
        through="Reactant",
        through_fields=("reaction", "compound"),
        related_name="_reactant_reactions",
    )

    _shorthand = "R"
    _name_field = "type"


class ProductModel(AbstractModel):
    class Meta:
        app_label = "hippo"
        abstract = True
        unique_together = ("reaction", "compound")

    amount = models.FloatField(default=1.0, validators=[MinValueValidator(0.0)])

    reaction = models.ForeignKey(
        "Reaction",
        on_delete=models.CASCADE,
        related_name="_products",
    )

    compound = models.ForeignKey(
        "Compound", on_delete=models.RESTRICT, related_name="_reaction_products"
    )

    solvent = models.ForeignKey(
        "Solvent",
        null=True,
        blank=True,
        on_delete=models.RESTRICT,
        related_name="_products",
    )


class ReactantModel(AbstractModel):
    class Meta:
        app_label = "hippo"
        abstract = True
        unique_together = ("reaction", "compound")

    amount = models.FloatField(default=1.0, validators=[MinValueValidator(0.0)])
    reaction = models.ForeignKey(
        "Reaction", on_delete=models.CASCADE, related_name="_reactants"
    )
    compound = models.ForeignKey(
        "Compound", on_delete=models.RESTRICT, related_name="_reaction_reactants"
    )

    solvent = models.ForeignKey(
        "Solvent",
        null=True,
        blank=True,
        on_delete=models.RESTRICT,
        related_name="_reactants",
    )


class SolventModel(AbstractModel):
    class Meta:
        app_label = "hippo"
        abstract = True

    name = models.CharField(max_length=30, unique=True)

    _name_field = "name"


### ANNOTATION


class TagModel(AbstractModel):
    class Meta:
        app_label = "hippo"
        abstract = True

    name = models.CharField(max_length=60, blank=False, unique=True)
    type = models.ForeignKey("TagType", on_delete=models.RESTRICT, related_name="_tags")

    _name_field = "name"


class TagTypeModel(AbstractModel):
    class Meta:
        abstract = True

    name = models.CharField(max_length=30, blank=False, unique=True)

    _name_field = "name"


class InspirationModel(AbstractModel):
    class Meta:
        abstract = True

    original = models.ForeignKey("Pose", on_delete=models.CASCADE, related_name="+")

    derivative = models.ForeignKey("Pose", on_delete=models.CASCADE, related_name="+")

    _score_types = models.ManyToManyField(
        "InspirationScoreType",
        through="InspirationScore",
        through_fields=("inspiration", "score_type"),
        related_name="_inspirations",
    )


class InspirationScoreTypeModel(AbstractModel):
    class Meta:
        abstract = True

    name = models.CharField(max_length=30)
    method = models.CharField(max_length=30)
    unit = models.CharField(max_length=30)
    description = models.CharField(max_length=120)


class InspirationScoreModel(AbstractModel):
    class Meta:
        abstract = True

    inspiration = models.ForeignKey(
        "Inspiration", on_delete=models.CASCADE, related_name="_scores"
    )
    score_type = models.ForeignKey(
        "InspirationScoreType", on_delete=models.CASCADE, related_name="_values"
    )

    value = models.FloatField()


class SubsiteModel(AbstractModel):
    class Meta:
        app_label = "hippo"
        abstract = True
        unique_together = ("target", "name")

    target = models.ForeignKey(
        "Target", on_delete=models.CASCADE, related_name="_subsites"
    )
    name = models.CharField(max_length=30)
    description = models.CharField(max_length=120, blank=True, null=True)
    metadata = models.JSONField(default=dict, blank=True)

    # _poses M2M field defined on

    _poses = models.ManyToManyField(
        "Pose",
        through="Observation",
        through_fields=("subsite", "pose"),
        related_name="_subsites",
    )

    _name_field = "name"


class ObservationModel(AbstractModel):
    class Meta:
        app_label = "hippo"
        abstract = True
        unique_together = ("subsite", "pose")

    subsite = models.ForeignKey(
        "Subsite", on_delete=models.RESTRICT, related_name="_observations"
    )
    pose = models.ForeignKey(
        "Pose", on_delete=models.CASCADE, related_name="_observations"
    )
    atom_ids = models.JSONField(default=dict)  # TODO: validation
    metadata = models.JSONField(default=dict, blank=True)

    _shorthand = "O"


class PlacementModel(AbstractModel):
    class Meta:
        app_label = "hippo"
        abstract = True
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


### FILE MANAGEMENT


class FileModel(AbstractModel):

    class Meta:
        abstract = True

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
        # pkl.tgz
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


### PROJECT MANAGEMENT


class CampaignModel(AbstractModel):

    class Meta:
        abstract = True

    name = models.CharField(max_length=30, unique=True)
    _targets = models.ManyToManyField("Target", related_name="_campaigns", blank=False)

    _name_field = "name"


class IterationModel(AbstractModel):

    class Meta:
        abstract = True

    ITERATION_STATUS = [
        ("P", "PENDING"),
        ("D", "DESIGN"),
        ("M", "MAKE"),
        ("T", "TEST"),
        ("A", "ANALYSE"),
        ("F", "FINISHED"),
        ("C", "CANCELLED"),
    ]

    number = models.PositiveSmallIntegerField(validators=[MinValueValidator(1)])

    campaign = models.ForeignKey(
        "Campaign",
        on_delete=models.CASCADE,
        related_name="_iterations",
    )

    status = models.CharField(max_length=1, choices=ITERATION_STATUS)

    _name_field = "long_name"


### RECIPES


class Route(AbstractModel):

    _reactions = models.ManyToManyField(
        "Reaction",
        through="RouteStep",
        through_fields=("route", "reaction"),
        related_name="_structures",
    )


class RouteStep(AbstractModel):

    route = models.ForeignKey("Route", on_delete=models.CASCADE, related_name="_steps")
    reaction = models.ForeignKey(
        "Reaction", on_delete=models.RESTRICT, related_name="_steps"
    )

    number = models.PositiveSmallIntegerField()

    multiplier = models.FloatField(default=1.0)

    is_root = models.BooleanField()
    is_leaf = models.BooleanField()


class RandomRecipeGenerator(AbstractModel):

    _suppliers = models.ManyToManyField(
        "Supplier",
        # through="RecipeComponent",
        # through_fields=("recipe", "route"),
        related_name="+",
    )

    budget = models.DecimalField(
        validators=[MinValueValidator(0.0)],
        blank=True,
        null=True,
        max_digits=8,
        decimal_places=2,
    )
    currency = models.CharField(max_length=3, blank=True, choices=CURRENCIES)

    origin = models.ForeignKey(
        "Recipe", on_delete=models.RESTRICT, related_name="_generators"
    )


class Recipe(AbstractModel):

    generator = models.ForeignKey(
        "RandomRecipeGenerator", on_delete=models.RESTRICT, related_name="_recipes"
    )

    _routes = models.ManyToManyField(
        "Route",
        through="RecipeComponent",
        through_fields=("recipe", "route"),
        related_name="_recipes",
    )


class RecipeComponent(AbstractModel):

    recipe = models.ForeignKey(
        "Recipe", on_delete=models.CASCADE, related_name="_recipe_components"
    )
    route = models.ForeignKey(
        "Route", on_delete=models.RESTRICT, related_name="_recipe_components"
    )

    multiplier = models.FloatField(default=1.0)


### SCORING


class RecipeScore(AbstractModel):

    recipe = models.ForeignKey(
        "Recipe", on_delete=models.CASCADE, related_name="_scores"
    )
    scorer = models.ForeignKey(
        "RecipeScorer", on_delete=models.RESTRICT, related_name="_scores"
    )

    score = models.FloatField(
        validators=[MinValueValidator(0.0), MaxValueValidator(1.0)]
    )


class RecipeScorer(AbstractModel):

    _recipes = models.ManyToManyField(
        "Recipe",
        through="RecipeScore",
        through_fields=("scorer", "recipe"),
        related_name="_scorers",
    )

    _attributes = models.ManyToManyField(
        "ScoringAttribute",
        through="ScoringAttributeWeight",
        through_fields=("scorer", "attribute"),
        related_name="_scorers",
    )


class ScoringAttribute(AbstractModel):

    name = models.CharField(max_length=30, unique=True)

    attribute = models.CharField(max_length=90)

    _recipes = models.ManyToManyField(
        "Recipe",
        through="AttributeValue",
        through_fields=("attribute", "recipe"),
        related_name="+",
    )


class AttributeValue(AbstractModel):

    attribute = models.ForeignKey(
        "ScoringAttribute", on_delete=models.RESTRICT, related_name="_values"
    )
    recipe = models.ForeignKey("Recipe", on_delete=models.CASCADE, related_name="+")

    value = models.FloatField()


class ScoringAttributeWeight(AbstractModel):

    attribute = models.ForeignKey(
        "ScoringAttribute", on_delete=models.CASCADE, related_name="+"
    )
    scorer = models.ForeignKey(
        "RecipeScorer", on_delete=models.RESTRICT, related_name="_attribute_weights"
    )

    is_inversed = models.BooleanField()

    weight = models.FloatField(
        validators=[MinValueValidator(0.0), MaxValueValidator(1.0)]
    )
