from ..models import SupplierModel
from .supplier_set import SupplierSet


class Supplier(SupplierModel):

    _objects = SupplierSet.as_manager()

    def __str__(self) -> str:
        return f'"{self.name}"'
