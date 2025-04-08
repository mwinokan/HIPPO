from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings

import mrich
from mrich import print
from pathlib import Path

from rdkit.Chem import PandasTools
from hippo.io.register import register_fragalysis_target


class Command(BaseCommand):
    help = "Import a Target and its Poses from fragalysis"

    def add_arguments(self, parser):
        parser.add_argument("target_name", type=str)
        parser.add_argument("-tas", "--target_access_string", type=str)
        parser.add_argument("-d", "--download_directory", type=str)
        parser.add_argument("-m", "--metadata_csv", type=str)
        parser.add_argument("-a", "--aligned_dir", type=str)
        parser.add_argument("-s", "--stack", type=str, default="production")
        parser.add_argument("-t", "--token", type=str)

    def handle(
        self,
        target_name: str,
        target_access_string: str | None = None,
        download_directory: "str | Path | None" = None,
        metadata_csv: "str | Path | None" = None,
        aligned_dir: "str | Path | None" = None,
        stack: str = "production",
        token: str | None = None,
        *args,
        **kwargs,
    ):

        mrich.h1("load_fragalysis")

        mrich.var("target_name", target_name)
        mrich.var("target_access_string", target_access_string)
        mrich.var("download_directory", download_directory)
        mrich.var("metadata_csv", metadata_csv)
        mrich.var("aligned_dir", aligned_dir)
        mrich.var("stack", stack)
        mrich.var("token", token)

        register_fragalysis_target(
            target_name=target_name,
            target_access_string=target_access_string,
            download_directory=download_directory,
            metadata_csv=metadata_csv,
            aligned_dir=aligned_dir,
            stack=stack,
            token=token,
        )
