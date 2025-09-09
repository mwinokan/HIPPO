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
import plotly.express as px

# import tempfile
import numpy as np

# import logging


class Command(BaseCommand):
    help = "Import scores a CSV for a given target"

    def add_arguments(self, parser):
        parser.add_argument("target_id", type=int)
        parser.add_argument("type_1", type=int)
        parser.add_argument("type_2", type=int)

    def handle(
        self,
        target_id: int,
        type_1: bool = False,
        type_2: bool = False,
        *args,
        **kwargs,
    ):

        target = Target.objects.get(pk=target_id)
        mrich.var("target", target)

        type_1 = ScoreType.objects.get(pk=type_1)
        mrich.var("type_1", type_1)
        type_2 = ScoreType.objects.get(pk=type_2)
        mrich.var("type_2", type_2)

        data = {}

        ### COMPOUND SCORES

        type_1_compoundscores = CompoundScore.objects.filter(type=type_1, target=target)

        mrich.var("#type_1_compoundscores", type_1_compoundscores.count())

        for cs in type_1_compoundscores:

            key = cs.compound_id

            data.setdefault(key, {type_1.id: [], type_2.id: []})

            data[key][type_1.id].append(cs.value)

        type_2_compoundscores = CompoundScore.objects.filter(type=type_2, target=target)

        mrich.var("#type_2_compoundscores", type_2_compoundscores.count())

        for cs in type_2_compoundscores:

            key = cs.compound_id

            data.setdefault(key, {type_1.id: [], type_2.id: []})

            data[key][type_2.id].append(cs.value)

        ### POSE SCORES

        type_1_posescores = PoseScore.objects.filter(
            pose__placement__structure__target=target, type=type_1
        )

        mrich.var("#type_1_posescores", type_1_posescores.count())

        for ps in type_1_posescores:

            compound_id = ps.pose.placement.compound_id

            data.setdefault(compound_id, {type_1.id: [], type_2.id: []})

            data[compound_id][type_1.id].append(ps.value)

        type_2_posescores = PoseScore.objects.filter(
            pose__placement__structure__target=target, type=type_2
        )

        mrich.var("#type_2_posescores", type_2_posescores.count())

        for ps in type_2_posescores:

            compound_id = ps.pose.placement.compound_id

            data.setdefault(compound_id, {type_1.id: [], type_2.id: []})

            data[compound_id][type_2.id].append(ps.value)

        plot_data = []

        for key, value in data.items():

            d = dict(
                compound_id=key,
            )

            values = value[type_1.id]

            if not values:
                continue

            d[f"{type_1} (avg)"] = np.mean(values)
            d[f"{type_1} (std)"] = np.std(values)

            values = value[type_2.id]

            if not values:
                continue

            d[f"{type_2} (avg)"] = np.mean(values)
            d[f"{type_2} (std)"] = np.std(values)

            plot_data.append(d)

        mrich.var("#data points", len(plot_data))

        if not plot_data:
            mrich.error("Empty plot data")
            return

        df = pd.DataFrame(plot_data)
        df = df.set_index("compound_id")
        mrich.writing(f"plot_T{target.id}_ST{type_1.id}_ST{type_2.id}.csv")
        df.to_csv(f"plot_T{target.id}_ST{type_1.id}_ST{type_2.id}.csv")

        file, _ = File.objects.get_or_create(
            path=f"plot_T{target.id}_ST{type_1.id}_ST{type_2.id}.csv",
            format_type=".csv",
            purpose="OUTPUT",
            content_type="META",
        )
        file.save()

        fig = px.scatter(
            df,
            x=f"{type_1} (avg)",
            y=f"{type_2} (avg)",
            error_x=f"{type_1} (std)",
            error_y=f"{type_2} (std)",
        )

        mrich.writing(f"plot_T{target.id}_ST{type_1.id}_ST{type_2.id}.html")
        fig.write_html(f"plot_T{target.id}_ST{type_1.id}_ST{type_2.id}.html")
        file, _ = File.objects.get_or_create(
            path=f"plot_T{target.id}_ST{type_1.id}_ST{type_2.id}.html",
            format_type=".html",
            purpose="OUTPUT",
            content_type="PLOT",
        )
        file.save()

        # fig.write_png(f"plot_T{target.id}_ST{type_1.id}_ST{type_2.id}.png")
        # file, _ = File.objects.get_or_create(path=f"plot_T{target.id}_ST{type_1.id}_ST{type_2.id}.png", format_type=".png", purpose="OUTPUT", content_type="PLOT")
        # file.save()
