import mrich

from hippo_django_orm.database import Database

db = Database(
    name="/Users/tfb64483/Software/HIPPO/tests_django_orm/test_django_models.sqlite"
)

# from hippo_django_orm.models import Compound
from hippo_django_orm.compound import Compound
from hippo_django_orm.target import Target
from hippo_django_orm.pose import Pose

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
    pass

pose = Pose.objects.first()
mrich.var("pose", pose)

all_compounds = Compound._objects.all()
mrich.var("all_compounds", all_compounds)

mrich.print(dir(Compound))
