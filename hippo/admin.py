from django.contrib import admin

# import models

from .protein import *
from .compound import *
from .pose import *
from .quoting import *
from .annotation import *
from .interactions import *
from .chemistry import *
from .files import *
from .projects import *

# setup custom wrappers

from .orm.setup import setup_models

setup_models()

# register models

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
admin.site.register(Structure)
admin.site.register(Placement)
admin.site.register(TagType)
admin.site.register(File)
admin.site.register(PoseScore)
admin.site.register(PoseScoreType)
admin.site.register(CompoundScore)
admin.site.register(CompoundScoreType)
admin.site.register(Campaign)
admin.site.register(Iteration)
