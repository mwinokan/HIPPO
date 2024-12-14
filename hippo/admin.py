from django.contrib import admin

# Register your models here.

from .target import Target
from .compound import Compound, CompoundScore, CompoundScoreType

# from .pose import Pose
from .quote import Quote
from .tag import Tag, TagType
from .subsite import Subsite
from .interaction import Interaction
from .feature import Feature
from .solvent import Solvent
from .observation import Observation
from .reaction import Reaction
from .reactant import Reactant
from .product import Product
from .supplier import Supplier
from .inspiration import Inspiration, InspirationScore, InspirationScoreType

from .models import (
    Campaign,
    Iteration,
)

from .file import File
from .pose import Pose, PoseScore, PoseScoreType

from .structure import Structure
from .placement import Placement

admin.site.register(Target)
admin.site.register(Compound)
admin.site.register(Pose)
admin.site.register(Quote)
admin.site.register(Tag)
admin.site.register(Subsite)
admin.site.register(Interaction)
admin.site.register(Feature)
admin.site.register(Solvent)
admin.site.register(Observation)
admin.site.register(Reaction)
admin.site.register(Reactant)
admin.site.register(Product)
admin.site.register(Supplier)
admin.site.register(Campaign)
admin.site.register(Iteration)
admin.site.register(Structure)
admin.site.register(Placement)
admin.site.register(TagType)
admin.site.register(File)
admin.site.register(PoseScore)
admin.site.register(PoseScoreType)
admin.site.register(CompoundScore)
admin.site.register(CompoundScoreType)
