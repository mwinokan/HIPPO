from ..models import PoseModel
from .pose_set import PoseSet


class Pose(PoseModel):

    _objects = PoseSet.as_manager()
    _parent_module = "pose"

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

        return value
