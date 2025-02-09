from .protein import Target, Structure
from .compound import Compound, CompoundScore, CompoundScoreType
from .pose import Pose, PoseScore, PoseScoreType
from .quoting import Quote, Supplier
from .annotation import (
    Tag,
    TagType,
    Subsite,
    Inspiration,
    InspirationScore,
    InspirationScoreType,
    Observation,
    Placement,
)
from .interactions import Interaction, Feature
from .chemistry import Solvent, Reaction, Reactant, Product
from .resources import File, Link
from .projects import Campaign, Iteration

MODELS = [
    Target,
    Compound,
    Pose,
    Quote,
    Tag,
    TagType,
    Link,
    Subsite,
    Interaction,
    Feature,
    Observation,
    Solvent,
    Reaction,
    Reactant,
    Product,
    Supplier,
    Structure,
    Placement,
    File,
    PoseScore,
    PoseScoreType,
    CompoundScore,
    CompoundScoreType,
    Inspiration,
    InspirationScore,
    InspirationScoreType,
    Campaign,
    Iteration,
]
