"""Custom django model fields

https://docs.djangoproject.com/en/5.1/howto/custom-model-fields/
"""

from django.db import models
from rdkit.Chem import Mol


class MolField(models.BinaryField):
    """Store an rdkit.Chem.Mol in a chemicalite.MOL field"""

    description = "an rdkit.Chem.Mol object"

    def db_type(self, connection):
        return "MOL"

    def from_db_value(self, value, expression, connection):
        if not value:
            return None

        # convert from chemicalite.MOL to Mol.ToBinary() format (in lieu of mol_to_binary_mol function)
        value = value.removeprefix(b"MOL\x00")

        return Mol(value)

    def to_python(self, value):
        if isinstance(value, Mol):
            return value

        if value is None:
            return value

        return Mol(value)

    def get_prep_value(self, value):
        if isinstance(value, bytes):
            return value

        binary = value.ToBinary()

        # convert from Mol.ToBinary() to chemicalite.MOL format (in lieu of mol_from_binary_mol function)
        binary = b"MOL\x00" + binary

        return binary
