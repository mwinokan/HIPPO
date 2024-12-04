from ..models import QuoteModel
from .quote_set import QuoteSet


class Quote(QuoteModel):

    _objects = QuoteSet.as_manager()
