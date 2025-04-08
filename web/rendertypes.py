from enum import Enum


class FieldRenderType(Enum):

    HIDDEN = 0
    TABLE = 1
    CARD = 2
    TOGGLE_CARD = 3


class ContentRenderType(Enum):

    TEXT_NORMAL = 1
    TEXT_MONOSPACE = 2
    INSTANCE_PILL = 3
    LIST_INSTANCE_PILLS = 4
    DICT_TABLE = 5
    MOL_2D_SVG = 6
    BOOL = 7
    OPTION = 7


DEFAULTS = {
    "<class 'django.db.models.fields.BigAutoField'>": dict(
        type=FieldRenderType.TABLE,
        content=ContentRenderType.TEXT_MONOSPACE,
        copyable=True,
    ),
    "<class 'django.db.models.fields.CharField'>": dict(
        type=FieldRenderType.TABLE,
        content=ContentRenderType.TEXT_NORMAL,
        copyable=True,
    ),
    "<class 'django.db.models.fields.reverse_related.ManyToManyRel'>": dict(
        type=FieldRenderType.TOGGLE_CARD, content=ContentRenderType.INSTANCE_PILL
    ),
    "<class 'django.db.models.fields.reverse_related.ManyToOneRel'>": dict(
        type=FieldRenderType.TOGGLE_CARD, content=ContentRenderType.INSTANCE_PILL
    ),
    "<class 'django.db.models.fields.related.ManyToManyField'>": dict(
        type=FieldRenderType.TOGGLE_CARD, content=ContentRenderType.INSTANCE_PILL
    ),
    "<class 'django.db.models.fields.BooleanField'>": dict(
        type=FieldRenderType.TABLE,
        content=ContentRenderType.TEXT_MONOSPACE,
    ),
    "<class 'django.db.models.fields.related.ForeignKey'>": dict(
        type=FieldRenderType.TABLE, content=ContentRenderType.INSTANCE_PILL
    ),
    "<class 'django.db.models.fields.reverse_related.OneToOneRel'>": dict(
        type=FieldRenderType.TABLE, content=ContentRenderType.INSTANCE_PILL
    ),
    "<class 'django.db.models.fields.TextField'>": dict(
        type=FieldRenderType.TOGGLE_CARD,
        content=ContentRenderType.TEXT_MONOSPACE,
        copyable=True,
    ),
}
