from django.db import models


class Compound(models.Model):

    compound_id = models.BigAutoField(primary_key=True)
    compound_inchikey = models.CharField(max_length=27)
    compound_alias = models.CharField(max_length=60)
    compound_smiles = models.CharField(max_length=90)
    # compound_base = models.ForeignKey(Compound, on_delete=models.CASCADE)
    compound_mol = models.BinaryField()
    compound_pattern_bfp.BinaryField()
    compound_morgan_bfp.BinaryField()
    compound_metadata.JSONField()
