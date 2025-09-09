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
import tempfile
import logging

mrich.debug("from fragmenstein import Wictor, Laboratory")
from fragmenstein import Wictor, Laboratory  # , Igor

mrich.debug("from fragmenstein.laboratory.validator import place_input_validator")
from fragmenstein.laboratory.validator import place_input_validator


class Command(BaseCommand):
    help = "Import scores a CSV for a given target"

    def add_arguments(self, parser):
        parser.add_argument("target_id", type=int)
        parser.add_argument("score_type_id", type=int)

    def handle(
        self,
        target_id: int,
        score_type_id: bool = False,
        n_retries: int = 3,
        n_cores: int = 8,
        timeout: int = 300,
        *args,
        **kwargs,
    ):

        target = Target.objects.get(pk=target_id)
        mrich.var("target", target)

        score_type = ScoreType.objects.get(pk=score_type_id)
        mrich.var("score_type", score_type)

        poses = (
            Pose.objects.filter(placement__structure__target=target)
            .exclude(scores__type=score_type)
            .distinct()
            .order_by("?")
        )

        n = poses.count()
        mrich.var("#poses", n)

        score_objs = []

        try:
            for i, pose in enumerate(poses):

                mrich.h1(f"{i} / {n}")

                with tempfile.TemporaryDirectory() as tmpdir:

                    scratch_dir = Path(tmpdir)

                    # set up lab
                    laboratory = setup_wictor_laboratory(
                        scratch_dir=scratch_dir,
                        pdb_block=pose.placement.structure.pdb_block,
                    )

                    # create inputs
                    queries = create_fragmenstein_queries_df(pose)

                    queries = place_input_validator(queries)

                    for attempt in range(n_retries):

                        # run the placement
                        result = laboratory.place(
                            queries,
                            n_cores=n_cores,
                            timeout=timeout,
                        )

                        # process outputs

                        if result is None:
                            mrich.error("Placement null result")
                            continue

                        assert len(result) == 1

                        result = result.iloc[0].to_dict()

                        if (
                            result["outcome"] == "crashed"
                            and result["error"] == "TimeoutError"
                        ):
                            mrich.error("Placement timed out")
                            continue

                        mrich.h3("Placement Result")

                        mrich.var("name", result.get("name", "N/A"))
                        mrich.var("error", result.get("error", "N/A"))
                        mrich.var("mode", result.get("mode", "N/A"))
                        mrich.var("∆∆G", result.get("∆∆G", "N/A"))
                        mrich.var("comRMSD", result.get("comRMSD", "N/A"))
                        mrich.var("runtime", result.get("runtime", "N/A"))
                        mrich.var("outcome", result.get("outcome", "N/A"))

                        break

                    if result and (value := result.get("∆∆G", None)):

                        d = dict(
                            pose=pose,
                            type=score_type,
                            value=float(value),
                            unit="kcal/mol",
                        )

                        if not value or pd.isna(value):
                            print(d)
                            continue

                        score_objs.append(PoseScore(**d))

        except KeyboardInterrupt:
            mrich.error("Interrupted!")
            pass

        mrich.bold(f"Inserting {len(score_objs)} PoseScore")
        PoseScore.objects.bulk_create(score_objs, ignore_conflicts=True)


def create_fragmenstein_queries_df(pose: "Pose"):

    return pd.DataFrame(
        [
            {
                "name": f"P{pose.id}",
                "smiles": pose.placement.compound.smiles,
                "hits": [pose.mol],
            }
        ]
    )


# from syndirella.slipper.SlipperFitter.setup_Fragmenstein
def setup_wictor_laboratory(
    *,
    scratch_dir: "Path",
    pdb_block: str,
    monster_joining_cutoff: float = 5,  # Å
) -> "Laboratory":

    # from fragmenstein import Laboratory, Wictor, Igor
    assert Path(scratch_dir).exists(), f"{scratch_dir=} does not exist"

    # set up Wictor
    Wictor.work_path = scratch_dir
    Wictor.monster_throw_on_discard = True  # stop if fragment unusable
    Wictor.monster_joining_cutoff = monster_joining_cutoff
    Wictor.quick_reanimation = False  # for the impatient
    Wictor.error_to_catch = Exception  # stop the whole laboratory otherwise
    Wictor.enable_stdout(logging.CRITICAL)
    Wictor.enable_logfile(scratch_dir / f"fragmenstein.log", logging.DEBUG)

    Laboratory.Victor = Wictor

    lab = Laboratory(pdbblock=pdb_block, covalent_resi=None, run_plip=False)

    return lab
