from django.apps import apps
from django.db.models import Manager, QuerySet
from rdkit import Chem

from .utils import registration_hash_tautomer_insensitive, superparent


class CompoundQueryset(QuerySet):
    def filter_qs(self):
        CompoundModel = apps.get_model('designdb', 'CompoundModel')
        qs = CompoundModel.objects.all()
        return qs

    def get_by_smiles(self, smiles):
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        try:
            sp = superparent(mol)
        except Exception as e:
            raise ValueError(f'SuperParent failed: {e}') from e

        h = registration_hash_tautomer_insensitive(sp)

        return self.filter_qs().get(compound_hash=h)


class CompoundManager(Manager):
    def get_queryset(self):
        return CompoundQueryset(self.model, using=self._db)

    # probably not needed here?
    def get_by_smiles(self, smiles):
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        try:
            sp = superparent(mol)
        except Exception as e:
            raise ValueError(f'SuperParent failed: {e}') from e

        h = registration_hash_tautomer_insensitive(sp)

        return self.get_queryset().get(compound_hash=h)
