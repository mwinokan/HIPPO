from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings

import mrich
from mrich import print
from pathlib import Path

from rdkit.Chem import PandasTools
from hippo.io.register import register_pose
from hippo.custom_models import Pose, Inspiration


class Command(BaseCommand):
    help = "Loads poses from an SDF"

    def add_arguments(self, parser):
        parser.add_argument("file", type=str)

    def handle(self, file: str, *args, **kwargs):

        mrich.h1("load_poses")

        mrich.var("file", file)

        file = Path(file)

        assert file.exists()

        match extension := file.name.split(".")[-1]:

            case "sdf":

                df = PandasTools.LoadSDF(str(file.resolve()))

            case _:

                raise NotImplementedError(f"unsupported file extension {extension}")

        mrich.var("#molecules", len(df))

        poses = {p.alias: p for p in Pose.objects.all()}
        inspiration_objects = []

        for i, row in df.iterrows():

            if row["ID"] == "ver_1.2":
                mrich.warning("Skipping Fragalysis header")
                continue

            if "ref_mols" in row:
                inspiration_names = row["ref_mols"].split(",")

                inspirations = []

                for inspiration_name in inspiration_names:

                    if inspiration_name in poses:
                        inspirations.append(poses[inspiration_name])
                    else:
                        s = f"Unknown inspiration: {inspiration_name}"
                        mrich.error(s)
                        raise ValueError(s)

            mol = row["ROMol"]

            data = dict(
                mol=mol,
            )

            pose = register_pose(**data)

            for inspiration in inspirations:
                inspiration_objects.append(
                    Inspiration(original=inspiration, derivative=pose)
                )

        Inspiration.objects.bulk_create(
            inspiration_objects,
            ignore_conflicts=True,
        )
