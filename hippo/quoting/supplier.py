from ..models import SupplierModel
from .supplier_set import SupplierSet


class Supplier(SupplierModel):

    _objects = SupplierSet.as_manager()
    _parent_module = "quoting"

    def __str__(self) -> str:
        return f'"{self.name}"'
