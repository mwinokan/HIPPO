"""This file defines the base models for the HIPPO Database schema using django, these classes are then subclassed in other files to create Python wrappers"""

from pathlib import Path
from datetime import date
from django.db import models
from .orm.expressions import MolFromSmiles, MolPatternBfpFromSmiles
from .orm.fields import MolField
from .orm.validators import validate_list_of_integers, validate_coord
from molparse.rdkit.features import FEATURE_FAMILIES, INTERACTION_TYPES
from django.core.validators import MinValueValidator, MaxValueValidator
from .abstract import AbstractModel, AbstractQuerySet


class TargetModel(AbstractModel):
    class Meta:
        app_label = "hippo"
        abstract = True

    name = models.CharField(max_length=60, blank=True, unique=True)
    metadata = models.JSONField(default=dict(), blank=True)

    _shorthand = "T"


class CompoundModel(AbstractModel):
    class Meta:
        app_label = "hippo"
        abstract = True

    inchikey = models.CharField(max_length=27, unique=True)
    alias = models.CharField(max_length=60, blank=True, unique=True, null=True)
    smiles = models.CharField(max_length=90)
    metadata = models.JSONField(default=dict(), blank=True)

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

    _scaffolds = models.ManyToManyField("Compound", related_name="_elaborations")
    _tags = models.ManyToManyField("Tag", related_name="_compounds")

    _shorthand = "C"


class PoseModel(AbstractModel):
    class Meta:
        app_label = "hippo"
        abstract = True

    inchikey = models.CharField(max_length=27, blank=True)
    alias = models.CharField(max_length=60, blank=True, unique=True, null=True)
    smiles = models.CharField(max_length=90, blank=True, null=True)

    reference = models.ForeignKey(
        "Pose", on_delete=models.SET_NULL, blank=True, null=True, related_name="+"
    )
    path = models.FilePathField(Path("/"), max_length=200, unique=True)
    compound = models.ForeignKey(
        "Compound", on_delete=models.CASCADE, related_name="_poses"
    )
    target = models.ForeignKey(
        "Target", on_delete=models.CASCADE, related_name="_poses"
    )
    mol = MolField(blank=True)
    energy_score = models.FloatField(blank=True, null=True)
    distance_score = models.FloatField(blank=True, null=True)
    metadata = models.JSONField(default=dict(), blank=True)
    fingerprinted = models.BooleanField(default=False)

    _inspirations = models.ManyToManyField("Pose", related_name="_derivatives")
    _tags = models.ManyToManyField("Tag", related_name="_poses")

    _shorthand = "P"


class QuoteModel(AbstractModel):
    class Meta:
        app_label = "hippo"
        abstract = True
        unique_together = ("amount", "supplier", "catalogue", "entry")

    from .price.price import CURRENCIES

    amount = models.FloatField(validators=[MinValueValidator(0.0)])
    supplier = models.CharField(max_length=60, blank=False)
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
        "Compound", on_delete=models.CASCADE, related_name="_quotes"
    )
    purity = models.FloatField(
        blank=True,
        null=True,
        validators=[MinValueValidator(0.0), MaxValueValidator(1.0)],
    )

    _shorthand = "Q"


class TagModel(AbstractModel):
    class Meta:
        app_label = "hippo"
        abstract = True

    name = models.CharField(max_length=60, blank=False, unique=True)

    _shorthand = None


class ReactionModel(AbstractModel):
    class Meta:
        app_label = "hippo"
        abstract = True

    type = models.CharField(max_length=60, blank=False)
    product = models.ForeignKey(
        "Compound", on_delete=models.CASCADE, related_name="_reactions"
    )
    yield_fraction = models.FloatField(
        validators=[MinValueValidator(0.0), MaxValueValidator(1.0)]
    )

    _shorthand = "R"


class ReactantModel(AbstractModel):
    class Meta:
        app_label = "hippo"
        abstract = True
        unique_together = ("reaction", "compound")

    amount = models.FloatField(validators=[MinValueValidator(0.0)])
    reaction = models.ForeignKey(
        "Reaction", on_delete=models.CASCADE, related_name="_reactants"
    )
    compound = models.ForeignKey(
        "Compound", on_delete=models.RESTRICT, related_name="+"
    )
    solvent = models.ForeignKey("Solvent", on_delete=models.SET_NULL, related_name="+")

    _shorthand = None


class FeatureModel(AbstractModel):

    class Meta:
        app_label = "hippo"
        abstract = True
        unique_together = (
            "family",
            "target",
            "chain_name",
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
    residue_number = models.PositiveSmallIntegerField()
    atom_names = models.JSONField()  # TODO: Validate as list

    _shorthand = "F"


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

    feature = models.ForeignKey(
        "Feature", on_delete=models.RESTRICT, related_name="_interactions"
    )
    atom_ids = models.JSONField(validators=[validate_list_of_integers])
    prot_coord = models.JSONField(validators=[validate_coord])
    lig_coord = models.JSONField(validators=[validate_coord])
    distance = models.FloatField(validators=[MinValueValidator(0.0)])
    angle = models.FloatField(
        blank=True,
        null=True,
        validators=[MinValueValidator(0.0), MaxValueValidator(360)],
    )
    energy = models.FloatField(blank=True, null=True)

    _shorthand = "I"


class SubsiteModel(AbstractModel):
    class Meta:
        app_label = "hippo"
        abstract = True
        unique_together = ("target", "name")

    target = models.ForeignKey(
        "Target", on_delete=models.CASCADE, related_name="_subsites"
    )
    name = models.CharField(max_length=30)
    metadata = models.JSONField(default=dict(), blank=True)
    _poses = models.ManyToManyField(
        "Pose",
        through="SubsiteMember",
        through_fields=("subsite", "pose"),
        related_name="_subsites",
    )

    _shorthand = "S"


class SubsiteMemberSet(AbstractQuerySet): ...


class SubsiteMember(AbstractModel):
    class Meta:
        app_label = "hippo"
        abstract = False
        unique_together = ("subsite", "pose")

    _objects = SubsiteMemberSet.as_manager()

    subsite = models.ForeignKey(
        "Subsite", on_delete=models.RESTRICT, related_name="_subsite_members"
    )
    pose = models.ForeignKey(
        "Pose", on_delete=models.CASCADE, related_name="_subsite_members"
    )
    atom_ids = models.JSONField()  # TODO: validation
    metadata = models.JSONField(default=dict(), blank=True)

    _shorthand = None


class SolventModel(AbstractModel):
    class Meta:
        app_label = "hippo"
        abstract = True

    name = models.CharField(max_length=30, unique=True)

    _shorthand = None
