from django.core.management.base import BaseCommand
from django.conf import settings
from django.utils import timezone
from django.db import transaction

import mrich
from mrich import print
from pathlib import Path

from hippo.tools import inchikey_from_smiles, sanitise_smiles
from hippo.models import *
import pandas as pd
import re


class Command(BaseCommand):
    help = "Import scores a CSV for a given target"

    def add_arguments(self, parser):
        parser.add_argument("upload_id", type=int)

    def handle(
        self,
        upload_id: int,
        debug: bool = False,
        *args,
        **kwargs,
    ):

        mrich.h1("Load Scores")
        mrich.var("upload_id", upload_id)

        upload = ScoreUpload.objects.get(id=upload_id)

        file_path = Path(upload.input_file.path)

        mrich.var("file_path", file_path)
        mrich.var("separator", upload.separator)
        mrich.var("identifier_column_name", upload.identifier_column_name)
        mrich.var("value_field", upload.value_field)
        mrich.var("unit_field", upload.unit_field)
        mrich.var("score_type", upload.score_type)
        mrich.var("target", upload.target)

        upload.status = 1
        upload.time_finished = None
        upload.message = "Registering File"

        file, _ = File.objects.get_or_create(
            path=file_path,
            format_type=".csv",
            content_type="META",
            purpose="INPUT",
        )

        upload.file = file
        upload.message = f"file: {file}"
        upload.save()

        try:

            ### LOAD DF

            df = pd.read_csv(file_path, sep=upload.separator)

            identifiers = []
            scores = []
            units = []

            match upload.related_model:
                case 0:

                    for i, row in df.iterrows():

                        identifiers.append(row[upload.identifier_column_name])

                        try:
                            scores.append(float(row[upload.value_field]))
                        except Exception as e:
                            mrich.warning(e)
                            continue

                        units.append(row[upload.unit_field])

                case _:
                    raise NotImplementedError("Pose score upload not yet supported")

            ### COMPOUNDS

            smiles_map = Compound.bulk_register(smiles=identifiers)

            upload.message += "\nCompounds registered"
            upload.save()

            compounds = {
                c.smiles: c
                for c in Compound.objects.filter(smiles__in=smiles_map.values())
            }
            compounds = {
                s_old: compounds[smiles] for s_old, smiles in smiles_map.items()
            }

            ### SCORES

            objs = []

            for s_old, score, unit in zip(identifiers, scores, units):

                compound = compounds[s_old]

                data = dict(
                    compound=compound,
                    type=upload.score_type,
                    target=upload.target,
                    value=score,
                    unit=unit,
                )

                objs.append(CompoundScore(**data))

            CompoundScore.objects.bulk_create(objs)

            upload.message += "\nScores registered"
            upload.save()

        #     with transaction.atomic():

        #         ### POSES

        #         poses = {}
        #         for i, row in df.iterrows():

        #             pose = Pose(
        #                 origin=upload.pose_origins,
        #                 smiles=row["smiles"],
        #                 inchikey=inchikey_from_smiles(row["smiles"]),
        #                 mol=row["ROMol"],
        #                 inspiration_distance=row[
        #                     upload.inspiration_distance_field_name
        #                 ],
        #                 binding_energy=row[upload.binding_energy_field_name],
        #                 ligand_energy=Pose.calculate_ligand_energy(row["ROMol"]),
        #             )

        #             poses[i] = pose

        #         ### PLACEMENTS

        #         placements = []

        #         for i, row in df.iterrows():

        #             structure_key = row[upload.protein_field_name]

        #             structure = structures[structure_key]

        #             assert (
        #                 structure
        #             ), f"No structure found for {upload.protein_field_name}={structure_key}"

        #             placements.append(
        #                 Placement(
        #                     structure=structure,
        #                     pose=poses[i],
        #                     compound=compounds[row["smiles"]],
        #                 )
        #             )

        #         Placement.objects.bulk_create(
        #             placements,
        #             ignore_conflicts=True,
        #             unique_fields=["structure", "pose"],
        #         )

        #         upload.message += "\nPlacements registered"
        #         upload.save()

        #         ### INSPIRATIONS

        #         pose_lookup = {pt.tag.name: pt.pose for pt in pose_query}

        #         inspirations = []

        #         for i, row in df.iterrows():

        #             ref_names = row[upload.inspirations_field_name].split(",")

        #             for ref_name in ref_names:

        #                 inspirations.append(
        #                     Inspiration(
        #                         original=pose_lookup[ref_name],
        #                         derivative=poses[i],
        #                     )
        #                 )

        #         Inspiration.objects.bulk_create(
        #             inspirations,
        #             ignore_conflicts=True,
        #             unique_fields=["original", "derivative"],
        #         )

        #         upload.message += "\nInspirations registered"
        #         upload.save()

        #         ### POSESETMEMBERS

        #         members = []
        #         for i, row in df.iterrows():

        #             members.append(
        #                 PoseSetMember(
        #                     parent=pset,
        #                     pose=poses[i],
        #                 )
        #             )

        #         PoseSetMember.objects.bulk_create(
        #             members, ignore_conflicts=True, unique_fields=["parent", "pose"]
        #         )

        #         upload.message += "\nPoseSetMembers registered"
        #         upload.save()

        #     ### UMAP

        #     if upload.compute_embedding:

        #         upload.message += "\nCalculating UMAP"
        #         upload.save()

        #         ok = pset.compute_embedding()

        #         if not ok:
        #             upload.message += "\nWARNING: Failed to compute UMAP"
        #             upload.save()

        #     upload.message += f"\nSUCCESS: Finished loading {upload}"
        #     upload.time_finished = timezone.now()
        #     upload.status = 2
        #     upload.save()

        except Exception as e:
            mrich.error(e)
            upload.message += f"\n{e.__class__.__name__.upper()}: {e}"
            upload.time_finished = timezone.now()
            upload.status = 3
            upload.save()
            raise
