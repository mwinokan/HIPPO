from django.core.management.base import BaseCommand
from django.conf import settings
from django.utils import timezone
from django.db import transaction

import mrich
from mrich import print
from pathlib import Path

from rdkit.Chem import PandasTools
from fragalysis.requests import download_target
from hippo.tools import inchikey_from_smiles
from hippo.models import *
import pandas as pd
import re


class Command(BaseCommand):
    help = "Import a Target and its Poses from fragalysis"

    def add_arguments(self, parser):
        parser.add_argument("upload_id", type=int)
        # parser.add_argument("--existing", type=str)
        # parser.add_argument("--debug", action="store_true")

    def handle(
        self,
        upload_id: int,
        # existing: str = None,
        # clear_previous_tag: bool = True,
        debug: bool = False,
        *args,
        **kwargs,
    ):

        mrich.h1("load_fragalysis")
        mrich.var("upload_id", upload_id)

        upload = SdfUpload.objects.get(id=upload_id)

        upload.status = 1
        upload.time_finished = None
        upload.message = "Registering File"
        upload.save()

        file_path = Path(upload.input_file.path)

        # register the File

        file, _ = File.objects.get_or_create(
            path=file_path,
            format_type=".sdf",
            content_type="LIGAND",
            purpose="INPUT",
        )

        # register the PoseSet

        upload.message += "\nRegistering PoseSet"
        upload.save()

        if not (pset := upload.pose_set):
            pset = PoseSet(name=file_path.name)
            pset.full_clean()
            pset.save()

            upload.pose_set = pset
            upload.save()

        ### LOAD DF

        df = PandasTools.LoadSDF(str(file_path.resolve()))

        upload.message += f"\nLoaded DataFrame, len={len(df)}"
        upload.save()

        ### PLACEMENTS

        ### COMPOUNDS

        ### POSES

        ### INSPIRATIONS

        ### POSESETMEMBERS

        print(df.iloc[0])

        # for i,row in df.iterrows():

        ### LOAD POSES FOR INSPIRATIONS

        # existing_poses = {
        #     p.alias: p
        #     for p in Pose.objects.prefetch_related(
        #         "placement__structure__target"
        #     ).filter(
        #         placement__structure__target=upload.target,
        #     )
        # }

        # pose_data = []
        # inspiration_objects = []

        # for i, row in df.iterrows():

        #     if row["ID"] == "ver_1.2":
        #         upload.message += "Skipping Fragalysis header"
        #         upload.save()
        #         continue

        #     mol = row["ROMol"]

        #     smiles =

        #     inchikey = inchikey_from_smiles(row["Smiles"])

        #     pose_data.append(dict(pose=Pose(

        #     )))

        # #     if "ref_mols" in row:
        # #         inspiration_names = row["ref_mols"].split(",")

        # #         inspirations = []

        # #         for inspiration_name in inspiration_names:

        # #             if inspiration_name in poses:
        # #                 inspirations.append(poses[inspiration_name])
        # #             else:
        # #                 s = f"Unknown inspiration: {inspiration_name}"
        # #                 mrich.error(s)
        # #                 raise ValueError(s)

        # #     mol = row["ROMol"]

        # #     data = dict(
        # #         mol=mol,
        # #     )

        # #     pose = register_pose(**data)

        # #     for inspiration in inspirations:
        # #         inspiration_objects.append(
        # #             Inspiration(original=inspiration, derivative=pose)
        # #         )

        # # Inspiration.objects.bulk_create(
        # #     inspiration_objects,
        # #     ignore_conflicts=True,
        # # )
