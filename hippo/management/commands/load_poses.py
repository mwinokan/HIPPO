from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings

import mrich
from mrich import print
from pathlib import Path

from rdkit.Chem import PandasTools
from hippo.io.register import register_pose


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

        database_settings = settings.DATABASES["default"]
        database_path = database_settings["NAME"]

        for i, row in df.iterrows():

            if row["ID"] == "ver_1.2":
                mrich.warning("Skipping Fragalysis header")
                continue

            mol = row["ROMol"]

            register_pose(mol=mol)
