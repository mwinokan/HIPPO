import mrich

from hippo_django_orm.database import Database
from rdkit.Chem import Mol
import molparse as mp

mrich.h3("TEST Database.__init__")

db = Database(
    name="/Users/tfb64483/Software/HIPPO/tests_django_orm/test_django_models.sqlite"
)

# from hippo_django_orm.models import Compound
from hippo_django_orm.compound import Compound
from hippo_django_orm.target import Target
from hippo_django_orm.pose import Pose

mrich.h3("TEST INSERT compound")

try:
    aspirin = Compound(
        inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        smiles="CC(=O)OC1=CC=CC=C1C(=O)O",
        alias="aspirin",
    )

    aspirin.full_clean()
    aspirin.save()
except Exception as e:
    mrich.error(e)
    pass

aspirin = Compound._objects.first()
assert aspirin
mrich.var("aspirin", aspirin)

mrich.h3("TEST INSERT target")

try:
    A71EV2A = Target(
        name="A71EV2A",
    )

    A71EV2A.full_clean()
    A71EV2A.save()
except Exception as e:
    mrich.error(e)
    pass

A71EV2A = Target.objects.first()
assert A71EV2A
mrich.var("A71EV2A", A71EV2A)

mrich.h3("TEST INSERT pose")

try:
    pose = Pose(
        # name="A71EV2A",
        compound=aspirin,
        target=A71EV2A,
        path="/path/to/pose.mol",
    )

    pose.full_clean()
    pose.save()
except Exception as e:
    mrich.error(e)
    raise
    pass

pose = Pose.objects.first()
mrich.var("pose", pose)

### TEST GETTING compound.mol

mrich.h3("TEST GETTING compound.mol")

compound = Compound._objects.first()
mrich.var("compound.mol", compound.mol)

### TEST GETTING pose.mol

mrich.h3("TEST GETTING pose.mol")

pose = Pose.objects.first()
mrich.var("pose.mol", pose.mol)
mrich.var("bool(pose.mol)", bool(pose.mol))

### TEST SETTING pose.mol

mrich.h3("TEST SETTING pose.mol")

mol = mp.rdkit.mol_from_smiles("CCC")
mrich.var("mol", mol)
pose.mol = mol
pose.full_clean()
pose.save()

### TEST GETTING pose.mol

mrich.h3("TEST GETTING pose.mol")

pose = Pose.objects.first()
mrich.var("pose.mol", pose.mol)
mrich.var("bool(pose.mol)", bool(pose.mol))

# pose.mol = mol

# pose.full_clean()
# pose.save()

# pose = Pose.objects.first()

# mrich.var("pose.mol", pose.mol)

# mrich.var("pdbblock", mp.rdkit.mol_to_pdb_block(pose.mol))

# all_compounds = Compound._objects.all()
# mrich.var("all_compounds", all_compounds)

# # mrich.print(dir(Compound))

# # mols = [m for m in Compound._objects.raw("SELECT id, mol_from_smiles(smiles) FROM compound")]

# # print(mols)

# # from django.db import connection

# # with connection.cursor() as cursor:
# #     cursor.execute("SELECT mol_to_binary_mol(mol_from_smiles(smiles)) FROM compound")
# #     mol, = cursor.fetchone()
# #     mrich.print(mol)
# #     mrich.print(Mol(mol))

# mrich.print(all_compounds[0].mol)
