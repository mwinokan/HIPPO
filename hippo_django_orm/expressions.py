"""Custom Expressions for GeneratedField

https://docs.djangoproject.com/en/5.1/ref/models/expressions/#django.db.models.Expression
"""

from django.db.models import Expression


class ChemicaliteFunction(Expression):
    """Create a Django Expression from any chemicalite function"""

    def __init__(self, expression, output_field):
        super().__init__(output_field=output_field)
        self.expression = expression

    def as_sql(self, compiler, connection, template=None):
        template = template or self.template
        data = {"expressions": self.expression}
        return template % data, []


class MolFromSmiles(ChemicaliteFunction):
    template = "mol_from_smiles( %(expressions)s )"


class MolPatternBfpFromSmiles(ChemicaliteFunction):
    template = "mol_pattern_bfp(mol_from_smiles( %(expressions)s ), 2048)"
