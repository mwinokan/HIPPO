import mrich
from rich import print
import hippo

legacy_animal = hippo.HIPPO("TestLegacyA71EV2A", "test_A71EV2A_legacy.sqlite")

# legacy_animal.add_hits(
#     target_name="A71EV2A",
#     metadata_csv="../tests/A71EV2A/metadata.csv",
#     aligned_directory="../tests/A71EV2A/aligned_files",
#     load_pose_mols=True,
# )

# from hippo_django_orm.database import Database
from hippo_django_orm import HIPPO

animal = HIPPO(
    "test_A71EV2A_django",
    "/Users/tfb64483/Software/HIPPO/tests_django_orm/test_A71EV2A_django.sqlite",
)

from hippo_django_orm.compound.compound import Compound
from hippo_django_orm.target.target import Target
from hippo_django_orm.pose.pose import Pose
from django.core.exceptions import ValidationError

### TARGETS

with mrich.loading("Target"):

    for legacy_target in legacy_animal.targets:

        target = Target(name=legacy_target.name)

        try:
            target.full_clean()
            target.save()
        except ValidationError as e:
            print(e)
            target = Target._objects.get(name=legacy_target.name)

### Compounds

with mrich.loading("Compound"):

    for legacy_compound in legacy_animal.compounds:
        # print(legacy_compound)
        # print(legacy_compound.dict)

        compound = Compound(
            smiles=legacy_compound.smiles,
            inchikey=legacy_compound.inchikey,
            alias=legacy_compound.alias,
            metadata=legacy_compound.metadata.data,
        )

        try:
            compound.full_clean()
            compound.save()
        except ValidationError as e:
            print(e)
            compound = Compound._objects.get(name=legacy_compound.name)

### Poses

with mrich.loading("Pose"):

    for legacy_pose in legacy_animal.poses:

        compound = Compound._objects.get(inchikey=legacy_compound.inchikey)

        # legacy_pose.summary()

        if legacy_pose.reference:
            reference = Pose._objects.get(path=legacy_pose.reference.path)
        else:
            reference = None

        inspirations = []

        if legacy_pose.inspirations:

            for legacy_inspiration in legacy_pose.inspirations:

                inspiration = Pose._objects.get(alias=legacy_inspiration.alias)
                mrich.print(inspiration)

                if inspiration:
                    inspirations.append(inspiration)

        target = Target._objects.get(name=legacy_pose.target.name)

        pose = Pose(
            alias=legacy_pose.alias,
            path=legacy_pose.path,
            compound=compound,
            target=target,
            reference=reference,
            # fingerprint=legacy_pose.has_fingerprint,
            energy_score=legacy_pose.energy_score,
            distance_score=legacy_pose.distance_score,
            metadata=legacy_pose.metadata.data,
        )  # , smiles=legacy_pose.smiles, mol=legacy_pose.mol)

        if inspirations:
            pose.inspirations.set(inspirations)

        try:
            pose.full_clean()
            pose.save()
        except ValidationError as e:
            print(e)
            pose = Pose._objects.get(alias=legacy_pose.alias)

        # break

### Test repr's

animal.summary()

mrich.var("compounds", animal.compounds)
mrich.var("poses", animal.poses)
mrich.var("targets", animal.targets)
mrich.var("targets", Target._objects.all())
mrich.var("compounds[:4]", Compound._objects.filter(id__lt=5))
mrich.var("len(compounds[:4])", len(Compound._objects.filter(id__lt=5)))

### Tags

## Fake Data

# Inspirations
# Scaffolds
# Reactions
