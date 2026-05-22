"""Main animal class for HIPPO"""

import logging
import re
from enum import Enum
from pathlib import Path

import mrich
import pandas as pd
from django.db import transaction

from .models import CompoundModel, PoseModel, TargetModel
from .services.ingestion import IngestionBatchResult, IngestionService
from .sets.compound import CompoundSet
from .sets.pose import PoseSet
from .utils import make_warn_once_per_key

logger = logging.getLogger(__name__)


class HIPPO:
    """Entry-point class of the xchem-hippo package.

    Update: this is atm not being called directly by the user.
    """

    def __init__(
        self,
        target_name: str,
    ) -> None:

        # TODO: user- or project based targets
        self._target, _ = TargetModel.objects.get_or_create(target_name=target_name)

        # TODO: the way this worked previously was it gave the HIPPO
        # instance full access to the pose table. When working with
        # multi-project central postgres db, this is almost certainly
        # not what I want. How is it that I'm going to keep this
        # updated? What does it mean upadte? Access to all objects
        # along this target?

        # self._compounds = CompoundTable(self.db)
        # self._poses = PoseSet(PoseModel.objects.all())  # <- NB! for testing
        # self._tags = TagTable(self.db)
        # self._reactions = ReactionTable(self.db)

        # ### in memory subsets
        # self._reactants = None
        # self._products = None
        # self._intermediates = None
        # self._scaffolds = None
        # self._elabs = None

    # @property
    # def name(self) -> str:
    #     """Returns the project name

    #     :returns: project name
    #     """
    #     return self._name

    @property
    def target(self) -> TargetModel:
        """Returns the target instance"""
        return self._target

    # actually expected to return all poses. filtering in PoseTable
    # class i.e. get_by_target.

    # Looks like I need to implement this. PoseService with some
    # manager- and instance mthods as helpers?

    # Actually it's more compplex than this: in the original code
    # there's PoseTable, and then there's PoseSet for a selection
    @property
    def poses(self):
        """Return pose instances for this target"""
        # return PoseModel.objects.filter(target=self._target)
        return PoseSet(PoseModel.objects.filter(target=self._target))

    @property
    def compounds(self) -> CompoundSet:
        """Return all compounds in the database"""
        return CompoundSet(CompoundModel.compound_filter.all())

    @property
    def num_poses(self) -> int:
        """Total number of Poses in the Database"""
        return self.poses.count()

    def add_hits(
        self,
        *,
        metadata_csv: str | Path,
        aligned_directory: str | Path,
        tags: list | None = None,
        skip: list | None = None,
        check_rmsd: bool = False,
        rmsd_threshold: float = 1.0,
        # debug: bool = False,
        # load_pose_mols: bool = False,
    ) -> pd.DataFrame:
        """Crystallographic hits from a Fragalysis download or XChemAlign alignment.

        For a Fragalysis download `aligned_directory` and `metadata_csv`
        should point to the `aligned_files` and `metadata.csv` at the
        root of the extracted download.
        For an XChemAlign dataset the `aligned_directory`
        should point to the `aligned_files`.

        :param target_name: Name of this protein :class:`.TargetModel`
        :param metadata_csv: Path to the metadata.csv from the Fragalysis download
        :param aligned_directory: Path to the aligned_files directory
            from the Fragalysis download
        :param skip: optional list of observation names to skip
        :param debug: bool:  (Default value = False)
        :returns: a DataFrame of metadata

        """

        ### Process arguments
        # NB! meta not required when loading XCA data
        assert metadata_csv, 'metadata.csv required'

        assert aligned_directory, 'aligned_directory must be provided'
        skip = skip or []
        tags = tags or ['hits']

        if not isinstance(aligned_directory, Path):
            aligned_directory = Path(aligned_directory)

        mrich.var('aligned_directory', aligned_directory)

        ### Determine data format

        # TODO: as it appears that users are currently only loading
        # fragalysis data, XCA format is not supported. Leaving the
        # format checks here to print a message for user

        class DataFormat(Enum):
            """DataFormat enum"""

            Fragalysis_v2 = 1
            XChemAlign_v2 = 2
            XChemAlign_v3 = 3

            def __str__(self) -> str:
                """name"""
                return self.name

        subdirs = list(aligned_directory.glob('*'))

        SUBDIR_PATTERN_FRAGALYSIS = re.compile(r'^.*\d{4}[a-z]$')
        SUBDIR_PATTERN_XCA = re.compile(r'^.*-.\d{4}$')

        fragalysis_subdirs_present = any(
            SUBDIR_PATTERN_FRAGALYSIS.match(subdir.name) for subdir in subdirs
        )
        xca_subdirs_present = any(
            SUBDIR_PATTERN_XCA.match(subdir.name) for subdir in subdirs
        )
        assert fragalysis_subdirs_present ^ xca_subdirs_present, (
            'Unexpected mixed data format'
        )

        if fragalysis_subdirs_present:
            data_format = DataFormat.Fragalysis_v2
        else:
            if any(list(subdir.glob('*_artefacts.pdb')) for subdir in subdirs):
                data_format = DataFormat.XChemAlign_v3
            else:
                data_format = DataFormat.XChemAlign_v2

            mrich.error(
                'Loading XChemAlign data currently not supported.'
                + ' Contact developers to enable this feature'
            )

        mrich.var('data_format', data_format)

        try:
            with transaction.atomic():
                result: IngestionBatchResult = IngestionService.ingest_filesystem(
                    root_path=aligned_directory,
                    target=self.target,
                    skip_records=skip,
                    compound_tag_list=tags,
                    metadata_file=metadata_csv,
                    check_rmsd=check_rmsd,
                    rmsd_threshold=rmsd_threshold,
                )
        except Exception as exc:
            logger.error(exc, exc_info=True)
            # TODO: handle gracefully
            raise Exception from exc

        # looking at the code, it seems to be the same, there are no
        # skips between observations and dirs_parsed declaratiosn
        mrich.var('#valid observations', result.attempts)

        # n_poses = self.num_poses
        # n_poses = PoseModel.objects.count()

        mrich.var('#directories parsed', result.attempts)
        mrich.var('#compounds registered', result.compounds_created)
        mrich.var('#poses registered', result.poses_created)

    def load_sdf(
        self,
        *,
        path: str | Path,
        reference: int | PoseModel | None = None,
        inspirations: list[int] | PoseSet | None = None,
        compound_tags: None | list[str] = None,
        pose_tags: None | list[str] = None,
        mol_col: str = 'ROMol',
        name_col: str = 'ID',
        inspiration_col: str = 'ref_mols',
        reference_col: str = 'ref_pdb',
        inspiration_map: None | dict = None,
        convert_floats: bool = True,
        skip_equal_dict: dict | None = None,
        skip_not_equal_dict: dict | None = None,
        check_rmsd: bool = False,
        rmsd_threshold: float = 1.0,
    ) -> None:
        """Add posed virtual hits from an SDF into the database.

        :param target: Name of the protein :class:`.TargetModel`
        :param path: Path to the SDF
        :param reference: Optional single reference :class:`.PoseModel` to use as
            the protein conformation for all poses, defaults to ``None``
        :param reference_col: Column that contains reference :class:`.PoseModel` aliases
            or ID's
        :param compound_tags: List of string Tags to assign to all created compounds,
            defaults to ``None``
        :param pose_tags: List of string Tags to assign to all created poses,
            defaults to ``None``
        :param mol_col: Name of the column containing the ``rdkit.ROMol`` ligands,
            defaults to ``"ROMol"``
        :param name_col: Name of the column containing the ligand name/alias,
            defaults to ``"ID"``
        :param inspirations: Optional single set of inspirations :class:`.PoseSet`
            object or list of IDs to assign as inspirations to all inserted poses,
            defaults to ``None``
        :param inspiration_col: Name of the column containing the list of inspiration
            :class:`.PoseModel` names or ID's, defaults to ``"ref_mols"``
        :param inspiration_map: Optional dictionary or callable mapping between
            inspiration strings found in ``inspiration_col`` and :class:`.PoseModel` ids
        :param energy_score_col: Name of the column containing the list of energy
            scores ``"energy_score"``
        :param distance_score_col: Name of the column containing the list of distance
            scores, defaults to ``"distance_score"``
        :param convert_floats: Try to convert all values to ``float``,
            defaults to ``True``
        :param skip_equal_dict: Skip rows where
            ``any(row[key] == value for key, value in skip_equal_dict.items())``,
            defaults to ``None``
        :param skip_not_equal_dict: Skip rows where
            ``any(row[key] != value for key, value in skip_not_equal_dict.items())``,
            defaults to ``None``

        All non-name columns are added to the PoseModel metadata.
        N.B. separate .mol files are not created. The molecule binary will only be
        stored in the .sqlite file and fake paths are added to the database.
        """
        # TODO: original code reads sdf into data frame. I don't see
        # much point for this in this function. get rid of it at some
        # point

        if not isinstance(path, Path):
            path = Path(path)

        skip_equal_dict = skip_equal_dict or {}
        skip_not_equal_dict = skip_not_equal_dict or {}

        mrich.debug(f'{path=}')

        compound_tags = compound_tags or []
        pose_tags = pose_tags or []

        if isinstance(inspirations, PoseSet):
            inspiration_list = list(inspirations.ids)
        elif isinstance(inspirations, list):
            # TODO: potentially check types
            inspiration_list = inspirations
        else:
            inspiration_list = []

        if reference and isinstance(reference, PoseModel):
            reference_id = reference.id
        else:
            reference_id = None

        if inspiration_map is None:
            inspiration_map = {}

        warn = make_warn_once_per_key()

        try:
            with transaction.atomic():
                result: IngestionBatchResult = IngestionService.ingest_sdf(
                    file_path=path,
                    target=self.target,
                    compound_tag_list=compound_tags,
                    pose_tag_list=pose_tags,
                    mol_col=mol_col,
                    name_col=name_col,
                    inspiration_col=inspiration_col,
                    inspirations=inspiration_list,
                    reference_col=reference_col,
                    reference=reference_id,
                    skip_equal=skip_equal_dict,
                    skip_not_equal=skip_not_equal_dict,
                    convert_floats=convert_floats,
                    field_warning=warn,
                    inspiration_map=inspiration_map,
                    check_rmsd=check_rmsd,
                    rmsd_threshold=rmsd_threshold,
                )
        except Exception as exc:
            logger.error(exc, exc_info=True)
            # TODO: handle gracefully
            raise Exception from exc

        # It's not clear what the original code was trying to do. I'm
        # going to issue warning when number of compounds and poses
        # was less than the number of compounds in sdf (not all were
        # successfully parsed) but that may not have been the original
        # intention
        if result.attempts == result.compounds_created:
            f = mrich.success
        else:
            f = mrich.warning

        f(f'{result.compounds_created} new compounds from {path}')

        if result.attempts == result.poses_created:
            f = mrich.success
        else:
            f = mrich.warning

        f(f'{result.poses_created} new poses from {path}')

    def add_syndirella_routes(
        self,
        pickle_path: str | Path,
        CAR_only: bool = True,
        pick_first: bool = True,
        check_chemistry: bool = True,
        register_routes: bool = True,
    ) -> pd.DataFrame:
        """Add routes found from syndirella --just_retro query"""

        try:
            with transaction.atomic():
                result: IngestionBatchResult = (
                    IngestionService.ingest_syndirella_routes(
                        pickle_path=pickle_path,
                        CAR_only=CAR_only,
                        pick_first=pick_first,
                        do_check_chemistry=check_chemistry,
                        register_routes=register_routes,
                    )
                )
        except Exception as exc:
            logger.error(exc, exc_info=True)
            # TODO: handle gracefully
            raise Exception from exc

    def add_syndirella_elabs(
        self,
        df_path: str | Path,
        max_energy_score: float | None = 0.0,
        max_distance_score: float | None = 2.0,
        require_intra_geometry_pass: bool = True,
        reject_flags: list[str] | None = None,
        register_reactions: bool = True,
        dry_run: bool = False,
        scaffold_route: 'RouteModel | None' = None,
        scaffold_compound: 'CompoundModel | None' = None,
        pose_tags: list[str] | None = None,
        product_tags: list[str] | None = None,
    ) -> pd.DataFrame:
        """
        Load Syndirella elaboration compounds and poses from a pickled DataFrame

        :param df_path: Path to the pickled DataFrame
        :param max_energy_score: Filter out poses with `∆∆G` above this value
        :param max_distance_score: Filter out poses with `comRMSD` above this value
        :param require_intra_geometry_pass: Filter out poses with falsy
            `intra_geometry_pass` values
        :param reject_flags: Filter out rows flagged with strings from this list
            (default = ["one_of_multiple_products",
            "selectivity_issue_contains_reaction_atoms_of_both_reactants"])
        :param scaffold_route: Supply a known single-step route to the scaffold product
            to use if scaffold placements are missing
        :param scaffold_compound: Supply a :class:`.CompoundModel` for the scaffold
            product to use if scaffold placements are missing
        :param dry_run: Don't insert new records into the database
            (for debugging/testing)
        :param pose_tags: Add these tags to all inserted poses, defaults to
            ["syndirella_product", "syndirella_placed"]
        :param product_tags: Add these tags to all inserted product compounds,
            defaults to ["syndirella_product"]
        :returns: annotated DataFrame
        """

        reject_flags = reject_flags or [
            'one_of_multiple_products',
            'selectivity_issue_contains_reaction_atoms_of_both_reactants',
        ]

        pose_tags = pose_tags or ['syndirella_product', 'syndirella_placed']
        product_tags = product_tags or ['syndirella_product']

        df_path = Path(df_path)
        mrich.h3(df_path.name)
        mrich.reading(df_path)
        df = pd.read_pickle(df_path)

        # testing
        # df = pd.read_csv(df_path.replace('.pkl.gz', '.csv'))

        try:
            with transaction.atomic():
                result: pd.DataFrame = IngestionService.ingest_syndirella_elabs(
                    df=df,
                    # TODO: check if target eists
                    target=self.target,
                    reject_flags=reject_flags,
                    pose_tag_list=pose_tags,
                    product_tag_list=pose_tags,
                    max_energy_score=max_energy_score,
                    max_distance_score=max_distance_score,
                    require_intra_geometry_pass=require_intra_geometry_pass,
                    register_reactions=register_reactions,
                    scaffold_route=scaffold_route,
                    scaffold_compound=scaffold_compound,
                )
                return result
        except Exception as exc:
            logger.error(exc, exc_info=True)
            # TODO: handle gracefully
            raise Exception from exc
