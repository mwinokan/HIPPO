from ..models import CompoundModel
from .compound_set import CompoundSet

# import mrich


class Compound(CompoundModel):
    class Meta:
        app_label = "hippo"

    _objects = CompoundSet.as_manager()
    _parent_module = "compound"

    def get_mol_svg_text(self, width=300, height=200):

        import re
        from rdkit.Chem.Draw import MolDraw2DSVG

        drawer = MolDraw2DSVG(width, height)
        drawer.DrawMolecule(self.mol)
        drawer.FinishDrawing()
        value = drawer.GetDrawingText()

        # transparent background
        value = re.sub(r"<rect style='opacity:1.0;fill:#FFFFFF.*> <\/rect>", "", value)

        return value
