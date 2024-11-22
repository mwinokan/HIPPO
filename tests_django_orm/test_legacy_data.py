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

from hippo_django_orm.compound import Compound
from hippo_django_orm.target import Target
from hippo_django_orm.pose import Pose
from hippo_django_orm.quote import Quote
from hippo_django_orm.interaction import Interaction
from hippo_django_orm.feature import Feature
from hippo_django_orm.subsite import Subsite
from hippo_django_orm.observation import Observation
from hippo_django_orm.reaction import Reaction
from hippo_django_orm.product import Product
from hippo_django_orm.reactant import Reactant

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

### Quote

quote = Quote(
    amount=0.1,
    supplier="Enamine",
    catalogue="EU BB",
    compound=animal.compounds[1],
)
quote.full_clean()
quote.save()

### Feature

t2 = Target(name="D68EV3C")
t2.full_clean()
t2.save()

feature = Feature(
    family="LumpedHydrophobe",
    target=animal.targets[1],
    chain_name="A",
    residue_name="CYS",
    residue_number="123",
    atom_names=["C", "N", "SG"],
)
feature.full_clean()
feature.save()

feature = Feature(
    family="LumpedHydrophobe",
    target=animal.targets[2],
    chain_name="A",
    residue_name="CYS",
    residue_number="123",
    atom_names=["C", "N", "SG"],
)
feature.full_clean()
feature.save()

### Interaction

features = animal.targets[1].features

mrich.print(type(features))
# mrich.print(dir(features))
mrich.print(list(f for f in features))

mrich.print(features.__class__)
mrich.print(animal.targets[1].features[0])

interaction = Interaction(
    family="LumpedHydrophobe",
    type="Hydrophobic",
    pose=animal.poses[1],
    feature=animal.targets[1].features[0],
    atom_ids=[1, 2, 3, 4],
    prot_coord=[1.0, 2.0, 3.0],
    lig_coord=[1.0, 2.0, 3.0],
    distance=0.2,
)
interaction.full_clean()
interaction.save()

### Subsite

subsite = Subsite(
    target=animal.T1,
    name="Oxyanion hole",
)


subsite.full_clean()
subsite.save()

observation = Observation(subsite=subsite, pose=animal.P1, atom_ids=[1, 2, 3])

observation.full_clean()
observation.save()

### Reaction

reaction = Reaction(
    type="Amidation",
)

reaction.full_clean()
reaction.save()

product = Product(reaction=reaction, compound=animal.C1)
product.full_clean()
product.save()

reactant = Reactant(reaction=reaction, compound=animal.C2)
reactant.full_clean()
reactant.save()

reactant = Reactant(reaction=reaction, compound=animal.C3)
reactant.full_clean()
reactant.save()

### Test repr's

animal.summary()

mrich.var("compounds", animal.compounds)
mrich.var("poses", animal.poses)
mrich.var("targets", animal.targets)
mrich.var("targets", Target._objects.all())
mrich.var("compounds[:4]", Compound._objects.filter(id__lt=5))
mrich.var("len(compounds[:4])", len(Compound._objects.filter(id__lt=5)))
mrich.var("targets[1].features", animal.targets[1].features)
mrich.var("animal.C1", animal.C1)
mrich.var("animal.T1", animal.T1)
mrich.var("animal.R1", animal.R1)
mrich.var("animal.T1.subsites", animal.T1.subsites)
mrich.var("T1.subsites[0].poses", animal.T1.subsites[0].poses)
mrich.var("T1.poses", animal.T1.poses, str(type(animal.T1.poses)))

### Tags

## Fake Data

# Inspirations
# Scaffolds
# Reactions
