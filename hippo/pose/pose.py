from ..models import PoseModel
from .pose_query_set import PoseQuerySet


class Pose(PoseModel):

    _objects = PoseQuerySet.as_manager()
    _parent_module = "pose"

    _custom_detail_view = True

    def get_mol_svg_text(self, width=300, height=200):

        import re
        from rdkit.Chem.Draw import MolDraw2DSVG
        from rdkit.Chem import MolToSmiles, MolFromSmiles

        mol = MolFromSmiles(MolToSmiles(self.mol))

        drawer = MolDraw2DSVG(width, height)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        value = drawer.GetDrawingText()

        # transparent background
        value = re.sub(r"<rect style='opacity:1.0;fill:#FFFFFF.*> <\/rect>", "", value)
        value = re.sub(r"<\?xml version='1\.0' encoding='iso-8859-1'\?>", "", value)

        return value
