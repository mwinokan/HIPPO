from django.core.management.base import BaseCommand
from django.conf import settings
from django.utils import timezone
from django.db import transaction

import mrich
from mrich import print
from pathlib import Path

from rdkit.Chem import PandasTools
from fragalysis.requests.download import download_target
from hippo.tools import inchikey_from_smiles
from hippo.models import *
import pandas as pd
import re
from numpy import nan


class Command(BaseCommand):
    help = "Import a Target and its Poses from fragalysis"

    def add_arguments(self, parser):
        parser.add_argument("download_id", type=int)
        parser.add_argument("--existing", type=str)
        parser.add_argument("--debug", action="store_true")

    def handle(
        self,
        download_id: int,
        existing: str = None,
        clear_previous_tag: bool = True,
        debug: bool = False,
        *args,
        **kwargs,
    ):

        mrich.h1("Load Fragalysis")
        mrich.var("download_id", download_id)
        mrich.var("existing", existing)

        destination = Path("data").resolve()

        download = FragalysisDownload.objects.get(id=download_id)

        try:
            download.status = 1
            download.message = "Requesting download"
            download.target = None
            download.time_finished = None
            download.save()

            if debug:
                download.summary()

            if existing:
                mrich.warning("Using existing download in", existing)
                target_dir = Path(existing)

            else:
                target_dir = download_target(
                    name=download.target_name,
                    tas=download.target_access_string,
                    stack=download.stack,
                    token=download.access_token,
                    destination=destination,
                )

            if not target_dir:
                download.status = 4
                download.message += "\nDownload failed"
                download.time_finished = timezone.now()
                download.save()
                return

            target_dir = target_dir.resolve()

            download.status = 2
            download.message += f"\n{target_dir=}"
            download.message += "\nDownload OK"
            download.save()

            # check for files

            if debug:
                mrich.debug("Checking files")

            try:
                metadata_csv = (target_dir / "metadata.csv").resolve()
                assert metadata_csv.exists(), "metadata.csv not in download"
                aligned_files = (target_dir / "aligned_files").resolve()
                assert aligned_files.exists(), "aligned_files not in download"
            except AssertionError as e:
                download.message += f"\n{e}"
                download.status = 4
                download.time_finished = timezone.now()
                return

            download.message += f"\n{metadata_csv=}"
            download.message += f"\n{aligned_files=}"
            download.message += "\nFiles OK"
            download.save()

            ### do the target loading

            ## TARGET

            if debug:
                mrich.debug("Creating Target")

            target, created = Target.objects.get_or_create(name=download.target_name)

            if created:
                download.message += f"\nCreated new target: {target}"
            else:
                download.message += f"\nUpdating existing target: {target}"
            download.save()

            df = pd.read_csv(metadata_csv)

            df["Compound code"] = df["Compound code"].replace({nan: None})

            missing = df[df["Smiles"] == "missing"]

            if len(missing):
                mrich.warning(len(df), "rows have 'missing' Smiles")
                df = df[df["Smiles"] != "missing"]

            ## COMPOUNDS

            if debug:
                mrich.debug("target", target)
                mrich.debug("Creating Compounds")

            smiles_map = Compound.bulk_register(smiles=df["Smiles"])

            download.message += "\nCompounds registered"

            compounds = {
                c.smiles: c
                for c in Compound.objects.filter(smiles__in=smiles_map.values())
            }

            compounds = {
                s_old: compounds[smiles]
                for s_old, smiles in smiles_map.items()
                if smiles in compounds
            }

            ## FILES

            if debug:
                mrich.debug("Creating Files")

            files = []

            for file in aligned_files.glob("*/*_apo-desolv.pdb"):
                files.append(
                    File(
                        path=file,
                        format_type=".pdb",
                        content_type="PROTEIN",
                        purpose="INPUT",
                    )
                )

            File.objects.bulk_create(
                files,
                ignore_conflicts=True,
                unique_fields=("path",),
            )

            files = {
                f.path: f
                for f in File.objects.filter(
                    format_type=".pdb",
                    content_type="PROTEIN",
                    purpose="INPUT",
                )
            }

            download.message += "\napo-desolv PDBs registered"

            ## STRUCTURES

            if debug:
                mrich.debug("Creating Structures")

            structures = []
            for short_code, alias in (
                df[["Code", "Experiment code"]]
                .drop_duplicates("Experiment code")
                .values
            ):

                p = str(aligned_files / short_code / (short_code + "_apo-desolv.pdb"))

                file = files[p]

                pdb_block = file.as_path.read_text(encoding="utf-8")

                structures.append(
                    Structure(
                        target=target,
                        alias=alias,
                        pdb_block=pdb_block,
                        protein_file=file,
                        origin="EXPERIMENT",
                        resolution=None,
                        metadata={},
                    )
                )

            Structure.objects.bulk_create(
                structures,
                update_conflicts=True,
                unique_fields=["protein_file"],
                update_fields=["pdb_block", "resolution", "metadata"],
            )

            structures = {
                s.alias: s
                for s in Structure.objects.filter(
                    target=target,
                    origin="EXPERIMENT",
                )
            }

            download.message += "\nStructures registered"

            with transaction.atomic():

                ## POSES

                if debug:
                    mrich.debug("Creating Poses")

                poses = {}

                for i, row in df.iterrows():

                    alias = row["Code"]
                    inchikey = inchikey_from_smiles(row["Smiles"])

                    sdf_file = aligned_files / alias / (alias + "_ligand.sdf")
                    sdf = PandasTools.LoadSDF(str(sdf_file))

                    if len(sdf) > 1:
                        download.message += (
                            f"\nWARNING: Multiple molecules in {sdf_file}"
                        )

                    mol = sdf["ROMol"].iloc[0]

                    pose = Pose(
                        inchikey=inchikey,
                        # alias=alias,
                        mol=mol,
                        smiles=row["Smiles"],
                        metadata={},
                        origin="EXPERIMENT",
                        ligand_energy=Pose.calculate_ligand_energy(mol),
                    )

                    poses[i] = pose

                Pose.objects.bulk_create(
                    poses.values(),
                    update_conflicts=True,
                    unique_fields=["mol"],
                    update_fields=["metadata"],
                )

                download.message += "\nPoses registered"
                download.save()

                ## PLACEMENTS

                if debug:
                    mrich.debug("Creating Placements")

                placements = []
                for i, row in df.iterrows():

                    structure = structures[row["Experiment code"]]
                    pose = poses[i]

                    compound = compounds.get(row["Smiles"], None)

                    if not compound:
                        download.message += f"\nERROR: Can't get compound with alias {row['Compound code']}"
                        download.message += f"\nERROR: Placement registration skipped ({structure}, {pose})"
                        continue

                    placements.append(
                        Placement(
                            structure=structure,
                            pose=pose,
                            compound=compound,
                        )
                    )

                Placement.objects.bulk_create(
                    placements,
                    ignore_conflicts=True,
                    unique_fields=["structure", "pose"],
                )

                download.message += "\nPlacements registered"
                download.save()

            ## PARSE TAGS

            # XCA_SITE_TAGS
            # XCA_SITE_TAG_FIELDS
            # CURATOR_TAG_TYPES

            pose_tags = []

            for col in df.columns:

                tag_name = None

                # observation code
                if col == "Code":
                    tag_type = "Observation Code"
                    tag_origin = "Fragalysis"

                # curator pose_tags
                elif match := re.search(r"^\[([a-zA-Z]*)\] (.*)$", col):
                    tag_type, tag_name = match.groups()
                    tag_origin = "Fragalysis Curator"

                elif match := re.search(r"(.*) alias", col):
                    (tag_type,) = match.groups()
                    tag_origin = "Fragalysis/XCA"

                else:
                    # mrich.warning(col)
                    continue

                for i, (pose_alias, value) in enumerate(df[["Code", col]].values):

                    if not value:
                        continue

                    if i not in poses:
                        continue

                    pose = poses[i]
                    pose_tags.append(
                        dict(
                            pose=pose,
                            tag_type=tag_type,
                            tag_name=tag_name or value,
                            tag_origin=tag_origin,
                        )
                    )

            comp_tags = []
            for comp_alias, smiles in df[["Compound code", "Smiles"]].values:

                if not comp_alias:
                    continue

                if comp_alias.startswith("Z"):
                    tag_type = "ZINC ID"
                elif comp_alias.startswith("ASAP-"):
                    tag_type = "ASAP ID"

                if smiles not in compounds:
                    continue

                comp_tags.append(
                    dict(
                        compound=compounds[smiles],
                        tag_type=tag_type,
                        tag_name=comp_alias,
                        tag_origin="Fragalysis/SoakDB",
                    )
                )

            pose_tag_df = pd.DataFrame(pose_tags)
            comp_tag_df = pd.DataFrame(comp_tags)
            combined_tag_df = pd.concat([pose_tag_df, comp_tag_df])

            ## CREATE TAG TYPES

            if debug:
                mrich.debug("Creating TagTypes")

            tag_types = {}
            for tag_type, tag_origin in (
                combined_tag_df[["tag_type", "tag_origin"]].drop_duplicates().values
            ):

                tt = TagType(
                    name=tag_type,
                    origin=tag_origin,
                )

                tag_types[(tag_type, tag_origin)] = tt

            TagType.objects.bulk_create(
                tag_types.values(),
                ignore_conflicts=True,
                unique_fields=["name", "origin"],
            )

            download.message += "\nTagTypes registered"
            download.save()

            ## CREATE TAGS

            if debug:
                mrich.debug("Creating Tags")

            tag_objs = []
            for tag_type, tag_origin, tag_name in (
                combined_tag_df[["tag_type", "tag_origin", "tag_name"]]
                .drop_duplicates()
                .values
            ):

                tag_type = tag_types[(tag_type, tag_origin)]

                if tag_type.id is None:
                    tag_type = TagType.objects.get(
                        name=tag_type.name, origin=tag_origin
                    )
                    tag_types[(tag_type.name, tag_origin)] = tag_type

                tag_objs.append(Tag(name=tag_name, type=tag_type))

            Tag.objects.bulk_create(
                tag_objs,
                ignore_conflicts=True,
                unique_fields=["name", "type"],
            )

            tag_objs = {(t.name, t.type_id): t for t in Tag.objects.all()}

            download.message += "\nTags registered"
            download.save()

            ## CREATE POSE TAGS

            if debug:
                mrich.debug("Creating PoseTags")

            delete_ids = set()
            pose_tag_objs = []
            for d in pose_tags:

                pose = d["pose"]
                tag_type = tag_types[d["tag_type"], d["tag_origin"]]
                tag_name = d["tag_name"]
                tag_obj = tag_objs[tag_name, tag_type.id]

                delete_ids.add(pose.id)

                pose_tag_objs.append(
                    PoseTag(
                        pose=pose,
                        tag=tag_obj,
                    )
                )

            with transaction.atomic():

                if clear_previous_tag:
                    PoseTag.objects.filter(pose_id__in=delete_ids).delete()

                PoseTag.objects.bulk_create(
                    pose_tag_objs,
                    ignore_conflicts=True,
                    unique_fields=["pose", "tag"],
                )

                download.message += "\nPoseTags registered"
                download.save()

            ## CREATE COMPOUND TAGS

            if debug:
                mrich.debug("Creating CompoundTags")

            delete_ids = set()
            # print(comp_tags)
            comp_tag_objs = []
            for d in comp_tags:

                compound = d["compound"]
                tag_type = tag_types[d["tag_type"], d["tag_origin"]]
                tag_name = d["tag_name"]
                tag_obj = tag_objs[tag_name, tag_type.id]

                delete_ids.add(compound.id)

                comp_tag_objs.append(
                    CompoundTag(
                        compound=compound,
                        tag=tag_obj,
                    )
                )

            with transaction.atomic():

                if clear_previous_tag:
                    CompoundTag.objects.filter(compound_id__in=delete_ids).delete()

                CompoundTag.objects.bulk_create(
                    comp_tag_objs,
                    ignore_conflicts=True,
                    unique_fields=["compound", "tag"],
                )
                download.message += "\nCompoundTags registered"
                download.save()

            ## SITES

            if debug:
                mrich.debug("Creating Sites")

            download.message += "\nLoading OK"
            download.status = 3
            download.time_finished = timezone.now()
            download.target = target
            download.save()

        except Exception as e:
            mrich.error(e)
            download.message += f"\n{e.__class__.__name__.upper()}: {e}"
            download.time_finished = timezone.now()
            download.status = 4
            download.save()
            raise
