from django.core.management.base import BaseCommand
from django.conf import settings
from django.utils import timezone
from django.db import transaction

import mrich
from mrich import print
from pathlib import Path

from rdkit.Chem import PandasTools, MolToSmiles
from hippo.tools import inchikey_from_smiles, sanitise_smiles
from hippo.models import *
import pandas as pd
import re


class Command(BaseCommand):
    help = "Import Poses from an SDF for a given target"

    def add_arguments(self, parser):
        parser.add_argument("upload_id", type=int)

    def handle(
        self,
        upload_id: int,
        debug: bool = False,
        *args,
        **kwargs,
    ):

        mrich.h1("Load SDF")
        mrich.var("upload_id", upload_id)

        upload = SdfUpload.objects.get(id=upload_id)

        upload.status = 1
        upload.time_finished = None
        upload.message = "Registering File"
        upload.save()

        file_path = Path(upload.input_file.path)

        try:

            ### LOAD DF

            df = PandasTools.LoadSDF(str(file_path.resolve()))

            if df.iloc[0]["ID"] == "ver_1.2":
                header = df.iloc[0].copy()
                df = df.iloc[1:].copy()
            else:
                header = None

            upload.message += f"\nLoaded DataFrame, len={len(df)}"
            upload.save()

            if upload.protein_field_name not in df.columns:
                upload.message += (
                    f"\nERROR: {upload.protein_field_name=} column required"
                )
                upload.status = 3
                upload.time_finished = timezone.now()
                upload.save()
                return

            if upload.inspirations_field_name not in df.columns:
                upload.message += (
                    f"\nERROR: {upload.inspirations_field_name=} column required"
                )
                upload.status = 3
                upload.time_finished = timezone.now()
                upload.save()
                return

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

            ### GET STRUCTURES

            structures = {
                s.alias: s for s in Structure.objects.filter(target=upload.target)
            }

            alias_tt = TagType.objects.get(name="Observation Code", origin="Fragalysis")

            pose_query = PoseTag.objects.select_related(
                "tag", "pose__placement__structure"
            ).filter(
                pose__placement__structure__target_id=upload.target_id,
                tag__type=alias_tt,
            )

            structures.update(
                {pt.tag.name: pt.pose.placement.structure for pt in pose_query}
            )

            ### COMPOUNDS

            df["smiles"] = df["ROMol"].apply(MolToSmiles)

            smiles_map = Compound.bulk_register(smiles=df["smiles"])

            upload.message += "\nCompounds registered"
            upload.save()

            compounds = {
                c.smiles: c
                for c in Compound.objects.filter(smiles__in=smiles_map.values())
            }
            compounds = {
                s_old: compounds[smiles] for s_old, smiles in smiles_map.items()
            }

            with transaction.atomic():

                ### POSES

                poses = {}
                for i, row in df.iterrows():

                    pose = Pose(
                        origin=upload.pose_origins,
                        smiles=row["smiles"],
                        inchikey=inchikey_from_smiles(row["smiles"]),
                        mol=row["ROMol"],
                    )

                    poses[i] = pose

                Pose.objects.bulk_create(
                    poses.values(),
                    update_conflicts=True,
                    unique_fields=["mol"],
                    update_fields=["origin"],
                )

                upload.message += "\nPoses registered"
                upload.save()

                ### PLACEMENTS

                placements = []

                for i, row in df.iterrows():

                    structure_key = row[upload.protein_field_name]

                    structure = structures[structure_key]

                    assert (
                        structure
                    ), f"No structure found for {upload.protein_field_name}={structure_key}"

                    placements.append(
                        Placement(
                            structure=structure,
                            pose=poses[i],
                            compound=compounds[row["smiles"]],
                        )
                    )

                Placement.objects.bulk_create(
                    placements,
                    ignore_conflicts=True,
                    unique_fields=["structure", "pose"],
                )

                upload.message += "\nPlacements registered"
                upload.save()

                ### INSPIRATIONS

                pose_lookup = {pt.tag.name: pt.pose for pt in pose_query}

                inspirations = []

                for i, row in df.iterrows():

                    ref_names = row[upload.inspirations_field_name].split(",")

                    for ref_name in ref_names:

                        inspirations.append(
                            Inspiration(
                                original=pose_lookup[ref_name],
                                derivative=poses[i],
                            )
                        )

                Inspiration.objects.bulk_create(
                    inspirations,
                    ignore_conflicts=True,
                    unique_fields=["original", "derivative"],
                )

                upload.message += "\nInspirations registered"
                upload.save()

                ### POSESETMEMBERS

                members = []
                for i, row in df.iterrows():

                    members.append(
                        PoseSetMember(
                            parent=pset,
                            pose=poses[i],
                        )
                    )

                PoseSetMember.objects.bulk_create(
                    members, ignore_conflicts=True, unique_fields=["parent", "pose"]
                )

                upload.message += "\nPoseSetMembers registered"
                upload.save()

            ### UMAP

            if upload.compute_umap:

                upload.message += "\nCalculating UMAP"
                upload.save()

                ok = pset.compute_umap()

                if not ok:
                    upload.message += "\nWARNING: Failed to compute UMAP"
                    upload.save()

            upload.message += f"\nSUCCESS: Finished loading {upload}"
            upload.time_finished = timezone.now()
            upload.status = 2
            upload.save()

        except Exception as e:
            # print(inchikey_from_smiles("Cn1nc(C(=O)NCCc2cocn2)c2ccccc21"))
            # print("Cn1nc(C(=O)NCCc2cocn2)c2ccccc21" in smiles_map.keys())
            # print("Cn1nc(C(=O)NCCc2cocn2)c2ccccc21" in compounds.keys())
            mrich.error(e)
            upload.message += f"\n{e.__class__.__name__.upper()}: {e}"
            upload.time_finished = timezone.now()
            upload.status = 3
            upload.save()
            raise
