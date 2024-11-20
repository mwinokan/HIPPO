from ..models import QuoteModel


class Quote(QuoteModel):

    def __str__(self):
        return f"Q{self.id}"
