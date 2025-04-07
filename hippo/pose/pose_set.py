from ..models import PoseSetModel
from .pose_sets import PoseSets
from .pose import Pose
from .pose_set_member import PoseSetMember
import mrich
from mrich import print


class PoseSet(PoseSetModel):

    _objects = PoseSets.as_manager()
    _parent_module = "pose"

    _custom_detail_view = True

    def calculate_pca(self):

        with mrich.loading("imports..."):
            from ..tools import get_cfps
            from sklearn.decomposition import PCA
            import numpy as np
            import pandas as pd

        members = list(self.poses.values("id", "pose"))

        pose_ids = [d["pose"] for d in members]

        poses = Pose.objects.filter(id__in=pose_ids)

        pose_mols = {d["id"]: d["mol"] for d in poses.values("id", "mol")}

        for member in members:
            mol = pose_mols[member["pose"]]
            member["mol"] = mol

        df = pd.DataFrame(members)

        with mrich.loading("Getting Compound fingerprints"):
            df["FP"] = df["mol"].map(get_cfps)

        X = np.array([x.fp for x in df["FP"]])

        with mrich.loading("Computing PCA"):
            pca = PCA(n_components=2, random_state=0)
            pca_fit = pca.fit_transform(X)

        df["PC1"] = pca_fit.T[0]
        df["PC2"] = pca_fit.T[1]

        members = []

        for i, row in df.iterrows():

            members.append(
                PoseSetMember(
                    id=row["id"],
                    pc1=row["PC1"],
                    pc2=row["PC2"],
                )
            )

        PoseSetMember.objects.bulk_update(members, fields=["pc1", "pc2"])

    def generate_tsnee_fig(self):

        import plotly.graph_objects as go

        n_valid = [v for v in self.poses.values_list("pc1", flat=True) if v is not None]

        if n_valid != self.poses.count():
            self.calculate_pca()

        members = self.poses

        fig = go.Figure()

        x = [m.pc1 for m in members]
        y = [m.pc2 for m in members]

        text = []
        colors = []
        pose_ids = []
        for m in members:

            pose = m.pose

            text.append(str(pose))
            pose_ids.append(pose.id)

            match m.review:
                case None:
                    colors.append("gray")
                case "GOOD":
                    colors.append("rgb(0,255,0)")
                case "BAD":
                    colors.append("red")

        trace = go.Scatter(
            x=x,
            y=y,
            mode="markers",
            marker=dict(
                size=12,
                color=colors,
            ),
            hovertemplate="%{text}<br>PC1: %{x}<br>PC2: %{y}<extra></extra>",
            text=text,
            customdata=pose_ids,
        )

        fig.add_trace(trace)

        fig.update_layout(
            margin=dict(l=20, r=20, t=20, b=20),
            paper_bgcolor="rgba(0,0,0,0)",
        )

        return fig
