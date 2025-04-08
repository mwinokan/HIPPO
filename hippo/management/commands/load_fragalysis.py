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


class Command(BaseCommand):
    help = "Import a Target and its Poses from fragalysis"

    def add_arguments(self, parser):
        parser.add_argument("download_id", type=int)

    def handle(
        self,
        download_id: int,
        *args,
        **kwargs,
    ):

        mrich.h1("load_fragalysis")

        destination = Path("data").resolve()

        download = FragalysisDownload.objects.get(id=download_id)
        download.summary()

        download.status = 1
        download.message = "Requesting download"
        download.target = None
        download.time_finished = None
        download.save()

        target_dir = download_target(
            name=download.target_name,
            tas=download.target_access_string,
            stack=download.stack,
            token=download.access_token,
            destination=destination,
        )

        # target_dir = Path("data/A71EV2A")

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

        target, created = Target.objects.get_or_create(name=download.target_name)

        if created:
            download.message += f"\nCreated new target: {target}"
        else:
            download.message += f"\nUpdating existing target: {target}"
        download.save()

        df = pd.read_csv(metadata_csv)

        ## COMPOUNDS

        Compound.bulk_register(smiles=df["Smiles"], alias=df["Compound code"])
        download.message += "\nCompounds registered"

        compounds = {c.alias: c for c in Compound.objects.all()}

        ## FILES

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

        structures = []
        for short_code, alias in (
            df[["Code", "Experiment code"]].drop_duplicates("Experiment code").values
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

            poses = []
            for i, row in df.iterrows():

                alias = row["Code"]
                inchikey = inchikey_from_smiles(row["Smiles"])

                sdf_file = aligned_files / alias / (alias + "_ligand.sdf")
                sdf = PandasTools.LoadSDF(str(sdf_file))

                if len(sdf) > 1:
                    download.message += f"\nWARNING: Multiple molecules in {sdf_file}"

                mol = sdf["ROMol"].iloc[0]

                poses.append(
                    Pose(
                        inchikey=inchikey,
                        alias=alias,
                        mol=mol,
                        smiles=row["Smiles"],
                        metadata={},
                        origin="EXPERIMENT",
                    )
                )

            Pose.objects.bulk_create(
                poses,
                update_conflicts=True,
                unique_fields=["alias", "mol"],
                update_fields=["metadata"],
            )

            poses = {
                p.alias: p
                for p in Pose.objects.filter(
                    origin="EXPERIMENT",
                )
            }

            ## PLACEMENTS

            placements = []
            for i, row in df.iterrows():

                structure = structures[row["Experiment code"]]
                pose = poses[row["Code"]]

                compound = compounds.get(row["Compound code"], None)

                if not compound:
                    download.message += (
                        f"\nERROR: Can't get compound with alias {row['Compound code']}"
                    )
                    download.message += f"\nERROR: Placement registration skipped ({structure=}, {pose=})"
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

        ## TAGS

        ## SITES

        download.message += "\nLoading OK"
        download.status = 3
        download.time_finished = timezone.now()
        download.target = target
        download.save()
