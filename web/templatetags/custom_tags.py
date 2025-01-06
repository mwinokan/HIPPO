from django import template

register = template.Library()


@register.filter("startswith")
def startswith(text, starts):
    if isinstance(text, str):
        return text.startswith(starts)
    return False


@register.filter("endswith")
def endswith(text, ends):
    if isinstance(text, str):
        return text.endswith(ends)
    return False


@register.filter("get")
def get(obj, attr):
    return getattr(obj, attr)
