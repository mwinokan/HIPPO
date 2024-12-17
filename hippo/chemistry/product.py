from ..models import ProductModel
from .product_set import ProductSet


class Product(ProductModel):

    _objects = ProductSet.as_manager()
    _parent_module = "chemistry"

    def __str__(self):
        return f"{self.reaction} -> Product_{self.id}"