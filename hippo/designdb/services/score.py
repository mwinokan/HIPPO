import logging
import re

# from mypackage.services.compound import CompoundService
# from rdkit.Chem import inchi
from designdb.models import PoseModel, ScoreValueModel, ScoringMethodModel

# from .validation.compound import ValidationError, validate_compound_data

SDF_XCAv2_PATTERN = re.compile(
    r'^.*-.\d{4}_._\d*_\d_.*-.\d{4}\+.\+\d*\+\d_ligand\.sdf$'
)
SDF_XCAV3_PATTERN = re.compile(
    r'^.*-.\d{4}_._\d*_._\d_.*-.\d{4}\+.\+\d*\+.\+\d_ligand\.sdf$'
)


SDF_FRAGALYSIS_PATTERN = re.compile(r'^.*\d{4}[a-z].sdf$')
PDBID_PATTERN = re.compile(r'^[A-Za-z0-9]{4}-[a-z].sdf$')


logger = logging.getLogger(__name__)


class ScoreService:
    def __init__(self, scoring_method_list: list[str] | None = None):
        self._scoring_method_list = scoring_method_list
        # self._score_map = {}
        self._scoring_method_cache = {}

        # unused, but I imagine this could take various arguments,
        # like include or exclude list

        # if self._scoring_method_list:
        #     for m in self._scoring_method_list:
        #         sm, _ = ScoringMethodModel.objects.get_or_create(
        #             method_name=m,
        #         )
        #         self._score_map[sm.method_name] = sm

    def add_scores_from_record(
        self,
        *,
        pose: PoseModel,
        record: dict[str, str | float],
        score_method_map: dict[str, ScoringMethodModel] | None = None,
    ):
        if score_method_map:
            scores = {col: record[col] for col in score_method_map if col in record}
        else:
            # FIXME: this because don't know how to select
            scores = {k: v for k, v in record.items() if k.lower().find('score') >= 0}

        for col_or_method_name, score_value in scores.items():
            if score_method_map:
                method = score_method_map[col_or_method_name]
            else:
                try:
                    method = self.scoring_methods[col_or_method_name]
                except KeyError:
                    method, _ = ScoringMethodModel.objects.get_or_create(
                        method_name=col_or_method_name,
                    )

            score = ScoreValueModel(
                pose=pose,
                compound=pose.compound,
                scoring_method=method,
                score={'score': score_value},
            )
            score.save()

    # def bulk_scores(poses: list[pose], record: dict[str, str | float]):
    #     # potentially lots of scores, can do bulk insertion all at once
    #     pass

    @property
    def scoring_methods(self) -> dict[str, ScoringMethodModel]:
        return self._scoring_method_cache
