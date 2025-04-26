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
        parser.add_argument("command", type=str)

    def handle(
        self,
        command,
        *args,
        **kwargs,
    ):

        mrich.var("command", command)
        exec(command)
