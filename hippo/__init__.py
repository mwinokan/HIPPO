"""

Hit Interaction Profiling for Progression Optimisation

HIPPO is a Python toolkit for structure- and fragment-based computational drug discovery,
storing large datasets in a database and facilitating rational decision making.
HIPPO was originally developed by Max Winokan while in XChem at Diamond Light Source.

See https://hippo-docs.winokan.com and https://github.com/mwinokan/HIPPO

"""

__version__ = "0.3.38"

from .animal import HIPPO
from .compound import Compound, Ingredient
from .cset import CompoundSet, CompoundTable, IngredientSet
from .db import Database
from .feature import Feature
from .metadata import MetaData
from .pose import Pose
from .price import Price
from .pset import PoseTable, PoseSet
from .quote import Quote
from .reaction import Reaction
from .recipe import Recipe, Route, RouteSet
from .rgen import RandomRecipeGenerator
from .rset import ReactionTable, ReactionSet
from .tags import TagTable, TagSet
from .target import Target
from .scoring import Scorer, CustomAttribute
from .web import ProjectPage
