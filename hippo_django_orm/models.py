"""This file defines the base models for the HIPPO Database schema using django, these classes are then subclassed in other files to create Python wrappers"""

from pathlib import Path
from django.db import models
from .expressions import MolFromSmiles, MolPatternBfpFromSmiles
from .fields import MolField
from datetime import date

from molparse.rdkit.features import FEATURE_FAMILIES, INTERACTION_TYPES


class TargetModel(models.Model):
    class Meta:
        app_label = "hippo"
        abstract = True

    id = models.BigAutoField(primary_key=True)
    name = models.CharField(max_length=60, blank=True, unique=True)
    metadata = models.JSONField(default=dict(), blank=True)


class CompoundModel(models.Model):
    class Meta:
        app_label = "hippo"
        abstract = True

    # regular fields

    id = models.BigAutoField(primary_key=True)
    inchikey = models.CharField(max_length=27, unique=True)
    alias = models.CharField(max_length=60, blank=True, unique=True, null=True)
    smiles = models.CharField(max_length=90)
    metadata = models.JSONField(default=dict(), blank=True)

    # chemicalite functions as custom Expressions in generated fields

    mol = models.GeneratedField(
        expression=MolFromSmiles("smiles", "mol"),
        output_field=MolField(blank=True),
        # output_field=models.BinaryField(blank=True),
        db_persist=True,
    )

    pattern_bfp = models.GeneratedField(
        expression=MolPatternBfpFromSmiles("smiles", "pattern_bfp"),
        output_field=models.BinaryField(blank=True),
        db_persist=True,
    )

    scaffolds = models.ManyToManyField("Compound", related_name="elaborations")
    tags = models.ManyToManyField("Tag", related_name="compounds")


class PoseModel(models.Model):
    class Meta:
        app_label = "hippo"
        abstract = True

    id = models.BigAutoField(primary_key=True)
    inchikey = models.CharField(max_length=27, blank=True)
    alias = models.CharField(max_length=60, blank=True, unique=True, null=True)
    smiles = models.CharField(max_length=90, blank=True, null=True)

    reference = models.ForeignKey(
        "Pose", on_delete=models.SET_NULL, blank=True, null=True, related_name="+"
    )
    path = models.FilePathField(Path("/"), max_length=200, unique=True)
    compound = models.ForeignKey(
        "Compound", on_delete=models.CASCADE, related_name="poses"
    )
    target = models.ForeignKey("Target", on_delete=models.CASCADE, related_name="poses")
    mol = MolField(blank=True)
    energy_score = models.FloatField(blank=True, null=True)
    distance_score = models.FloatField(blank=True, null=True)
    metadata = models.JSONField(default=dict(), blank=True)
    fingerprinted = models.BooleanField(default=False)

    inspirations = models.ManyToManyField("Pose", related_name="derivatives")
    tags = models.ManyToManyField("Tag", related_name="poses")


class QuoteModel(models.Model):
    class Meta:
        app_label = "hippo"
        abstract = True
        unique_together = ("amount", "supplier", "catalogue", "entry")

    from .price import CURRENCIES

    id = models.BigAutoField(primary_key=True)
    amount = models.FloatField()  # TODO: ADD Positive validation
    supplier = models.CharField(max_length=60, blank=False)
    catalogue = models.CharField(max_length=60, blank=True)
    entry = models.CharField(max_length=60, blank=True)
    lead_time = models.FloatField()  # TODO: ADD positive validation
    price = models.DecimalField(
        max_digits=8, decimal_places=2
    )  # TODO: ADD positive validation
    currency = models.CharField(max_length=3, blank=True, choices=CURRENCIES)
    date = models.DateField(default=date.today)
    compound = models.ForeignKey(
        "Compound", on_delete=models.CASCADE, related_name="quotes"
    )
    purity = models.FloatField()  # TODO: ADD 0-1 range validation


class TagModel(models.Model):
    class Meta:
        app_label = "hippo"
        abstract = True

    id = models.BigAutoField(primary_key=True)
    name = models.CharField(max_length=60, blank=False, unique=True)


class ReactionModel(models.Model):
    class Meta:
        app_label = "hippo"
        abstract = True

    id = models.BigAutoField(primary_key=True)
    type = models.CharField(max_length=60, blank=False)
    product = models.ForeignKey(
        "Compound", on_delete=models.CASCADE, related_name="reactions"
    )
    yield_fraction = models.FloatField()  # TODO: ADD validation


class ReactantModel(models.Model):
    class Meta:
        app_label = "hippo"
        abstract = True
        unique_together = ("reaction", "compound")

    id = models.BigAutoField(primary_key=True)
    amount = models.FloatField()  # TODO: validate positive
    reaction = models.ForeignKey(
        "Reaction", on_delete=models.CASCADE, related_name="reactants"
    )
    compound = models.ForeignKey("Compund", on_delete=models.RESTRICT, related_name="+")
    solvent = models.ForeignKey("Solvent", on_delete=models.SET_NULL, related_name="+")


class FeatureModel(models.Model):
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

    id = models.BigAutoField(primary_key=True)

    family = models.CharField(max_length=30, choices=FEATURE_FAMILIES)

    target = models.ForeignKey(
        "Target", on_delete=models.CASCADE, related_name="features"
    )

    chain_name = models.CharField(max_length=5)
    residue_name = models.CharField(max_length=10)
    residue_number = models.PositiveSmallIntegerField()
    atom_names = models.JSONField()  # TODO: Validate as list


class InteractionModel(models.Model):
    class Meta:
        app_label = "hippo"
        abstract = True
        unique_together = ("feature", "pose", "type", "family")

    id = models.BigAutoField(primary_key=True)

    family = models.CharField(max_length=30, choices=FEATURE_FAMILIES)
    type = models.CharField(max_length=30, choices=INTERACTION_TYPES.values())

    pose = models.ForeignKey(
        "Pose", on_delete=models.CASCADE, related_name="interactions"
    )

    feature = models.ForeignKey(
        "Feature", on_delete=models.RESTRICT, related_name="interactions"
    )
    atom_ids = models.JSONField()  # TODO: Validation
    prot_coord = models.JSONField()  # TODO: Validation
    lig_coord = models.JSONField()  # TODO: Validation
    distance = models.FloatField()  # TODO: Validate positive
    angle = models.FloatField()  # TODO: Validate positive
    energy = models.FloatField()


class SubsiteModel(models.Model):
    class Meta:
        app_label = "hippo"
        abstract = True
        unique_together = ("target", "name")

    id = models.BigAutoField(primary_key=True)
    target = models.ForeignKey(
        "Target", on_delete=models.CASCADE, related_name="subsites"
    )
    name = models.CharField(max_length=30)
    metadata = models.JSONField(default=dict(), blank=True)
    poses = models.ManyToManyField(
        "Pose",
        through="SubsiteMember",
        through_fields=("subsite", "pose"),
        related_name="subsites",
    )


class SubsiteMember(models.Model):
    class Meta:
        app_label = "hippo"
        abstract = False
        unique_together = ("subsite", "pose")

    id = models.BigAutoField(primary_key=True)
    subsite = models.ForeignKey(
        "Subsite", on_delete=models.RESTRICT, related_name="_subsite_members"
    )
    pose = models.ForeignKey(
        "Pose", on_delete=models.CASCADE, related_name="_subsite_members"
    )
    atom_ids = models.JSONField()  # TODO: validation
    metadata = models.JSONField(default=dict(), blank=True)


class SolventModel(models.Model):
    class Meta:
        app_label = "hippo"
        abstract = True

    id = models.BigAutoField(primary_key=True)
    name = models.CharField(max_length=30, unique=True)
