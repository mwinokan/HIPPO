import mrich
from django.core.exceptions import ValidationError
from django.utils.translation import gettext_lazy


def validate_list_of_integers(value):
    for v in value:
        try:
            int(v)
        except Exception as e:
            mrich.error(e)
            raise ValidationError(
                gettext_lazy("%(value)s is not a list of integers"),
                params={"value": value},
            )


def validate_coord(value):
    try:
        assert len(value) == 3
        for v in value:
            float(v)

    except Exception as e:
        mrich.error(e)
        raise ValidationError(
            gettext_lazy("%(value)s is not a 3D coordinate"),
            params={"value": value},
        )
