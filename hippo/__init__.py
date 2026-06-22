from .bootstrap import load_hippo as HIPPO

__all__ = ['HIPPO', 'IngredientSet']


def __getattr__(name):
    """Lazily expose select designdb classes at the package top level.

    Imported on first access (after ``load_hippo()``/``HIPPO()`` has configured
    Django) rather than at ``import hippo`` time, which runs before Django is
    configured. NB: transitional -- exposing internal classes like this is
    pending the client-exposure design (see RecipeManager discussion).
    """
    if name == 'IngredientSet':
        from designdb.sets.ingredient import IngredientSet

        return IngredientSet
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')
