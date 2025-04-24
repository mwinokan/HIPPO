from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
import mrich
from mrich import print
import pandas as pd

# from pathlib import Path

from hippo.models import *


class Command(BaseCommand):
    help = "Export human reviews of a given poseset"

    def add_arguments(self, parser):
        parser.add_argument("poseset_id", type=str)

    def handle(self, poseset_id, *args, **kwargs):

        mrich.h1("export_poseset_reviews")
        mrich.var("poseset_id", poseset_id)

        pset = PoseSet.objects.get(id=poseset_id)

        pset.summary()

        df = pd.DataFrame(pset.poses.values("pose_id", "pc1", "pc2", "review"))
        df.set_index("pose_id", inplace=True)
        df.to_csv(f"PoseSet_{poseset_id}.csv")
