"""This file defines the base models for the HIPPO Database schema using django, these classes are then subclassed in other files to create Python wrappers"""

from pathlib import Path
from django.db import models
from .expressions import MolFromSmiles, MolPatternBfpFromSmiles
from .fields import MolField


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
    alias = models.CharField(max_length=60, blank=True, unique=True)
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


class PoseModel(models.Model):
    class Meta:
        app_label = "hippo"
        abstract = True

    id = models.BigAutoField(primary_key=True)
    inchikey = models.CharField(max_length=27, unique=True, blank=True)
    alias = models.CharField(max_length=60, blank=True, unique=True)
    smiles = models.CharField(max_length=90, blank=True, unique=True)

    reference = models.ForeignKey(
        "Pose", on_delete=models.SET_NULL, blank=True, null=True
    )
    path = models.FilePathField(Path("/"), max_length=200)
    compound = models.ForeignKey("Compound", on_delete=models.CASCADE)
    target = models.ForeignKey("Target", on_delete=models.CASCADE)
    # mol = MolField(blank=True)
    mol = models.BinaryField(blank=True)
    energy_score = models.FloatField(blank=True, null=True)
    distance_score = models.FloatField(blank=True, null=True)
    metadata = models.JSONField(default=dict(), blank=True)

    ### legacy (pre-django implementation)
    # fingerprint
