"""This file defines the base models for the HIPPO Database schema using django, these classes are then subclassed in other files to create Python wrappers"""

from pathlib import Path
from django.db import models


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

    id = models.BigAutoField(primary_key=True)
    inchikey = models.CharField(max_length=27, unique=True)
    alias = models.CharField(max_length=60, blank=True, unique=True)
    smiles = models.CharField(max_length=90)
    mol = models.BinaryField(blank=True)
    metadata = models.JSONField(default=dict(), blank=True)

    ### legacy (pre-django implementation)
    # base = models.ForeignKey(Compound, on_delete=models.CASCADE)
    # pattern_bfp = models.BinaryField()
    # morgan_bfp = models.BinaryField()


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
    mol = models.BinaryField(blank=True)
    energy_score = models.FloatField(blank=True, null=True)
    distance_score = models.FloatField(blank=True, null=True)
    metadata = models.JSONField(default=dict(), blank=True)

    ### legacy (pre-django implementation)
    # fingerprint
