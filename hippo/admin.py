from django.contrib import admin

# Register your models here.

from .target import Target
from .compound import Compound
from .pose import Pose
from .quote import Quote
from .tag import Tag
from .subsite import Subsite
from .interaction import Interaction
from .feature import Feature
from .solvent import Solvent
from .observation import Observation
from .reaction import Reaction
from .reactant import Reactant
from .product import Product
from .supplier import Supplier

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
