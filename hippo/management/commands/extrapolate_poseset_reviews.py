from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command

import pandas as pd
import mrich
from mrich import print
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C

from hippo.models import *


class Command(BaseCommand):
    help = "extrapolate human reviews of a given poseset"

    def add_arguments(self, parser):
        parser.add_argument("poseset_id", type=str)

    def handle(self, poseset_id, *args, **kwargs):

        pset = PoseSet.objects.get(id=poseset_id)

        pset.extrapolate_reviews()
