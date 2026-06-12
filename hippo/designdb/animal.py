"""Main animal class for HIPPO"""

import logging
import re
from datetime import datetime
from enum import Enum
from pathlib import Path

import mrich
import pandas as pd
from django.db import transaction

from .models import (
    CompoundModel,
    EnumerationMethodModel,
    PoseMethodModel,
    PoseModel,
    Project,
    ScoringMethodModel,
    TargetModel,
)
from .services.download import DownloadService
from .services.ingestion import IngestionBatchResult, IngestionService
from .services.method import MethodService
from .services.route import RouteService
from .services.subsite import SubsiteService
from .sets.compound import CompoundSet
from .sets.pose import PoseSet
from .settings import DEFAULT_POSE_METHODS
from .utils import make_warn_once_per_key

logger = logging.getLogger(__name__)

# Root directory under which Fragalysis downloads are extracted, laid out as
# data/downloads/<project_name>/<target_name>.
DOWNLOADS_DIR = Path('data') / 'downloads'

# Fragalysis download flags requested when fetching the full hit data for a
# target (see HIPPO._ensure_hit_data). Everything else stays False.
HIT_DATA_FLAGS = (
    'apo_file',
    'bound_file',
    'apo_solv_file',
    'apo_desolv_file',
    'ligand_pdb',
    'ligand_sdf',
    'ligand_smiles',
    'sdf_info',
    'smiles_info',
    'metadata_info',
)


class HIPPO:
    """Entry-point class of the xchem-hippo package.

    Update: this is atm not being called directly by the user.
    """

    def __init__(
        self,
        target_name: str,
        target_access_string: str,
    ) -> None:

        # TODO: with working db, hippo shouldn't be creating projects
        project, _ = Project.objects.get_or_create(
            project_name=target_access_string,
        )

        # TODO: user- or project based targets
        self._target, _ = TargetModel.objects.get_or_create(
            target_name=target_name,
            project=project,
        )

        # Download state (see _ensure_hit_data / _ensure_apo_desolv_files). The
        # full hit data persists on disk and is reused across sessions; the
        # apo_desolv subset is re-downloaded once per instance to stay fresh.
        self._hit_data_path: Path | None = None
        self._apo_desolv_path: Path | None = None
        self._apo_desolv_downloaded_at: datetime | None = None

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

    def _ensure_hit_data(
        self, auth_token: str | None = None, stack: str = 'production'
    ) -> Path:
        """Ensure this target's full crystallographic hit data is available locally.

        Downloads the target's Fragalysis data (all observations) from the stack
        via :class:`.DownloadService` and returns the path to the extracted
        directory (``data/downloads/<project>/<target>``). The requested file
        types are :data:`.HIT_DATA_FLAGS` (apo/bound/ligand/sdf/smiles/metadata).

        Unlike :meth:`._ensure_apo_desolv_files`, this data **persists**: if a
        previous full download is already on disk it is reused without
        re-fetching (detected by the presence of ``metadata.csv``, which an
        apo_desolv-only download does not produce).

        .. note::
            HIPPO-level helper, intended to be called from user-facing
            :class:`.HIPPO` methods (e.g. :meth:`.add_hits`). It must not be
            called from the components or services layer.

        :param auth_token: optional Fragalysis ``sessionid``; otherwise the
            ``FRAGALYSIS_AUTH_TOKEN`` environment variable is used
        :param stack: Fragalysis stack to download from, a key into
            :data:`.STACK_URLS`; defaults to ``'production'``
        :returns: path to the extracted download directory
        """

        target_name = self._target.target_name
        project_name = self._target.project.project_name

        destination = DOWNLOADS_DIR / project_name
        target_dir = destination / target_name

        # reuse an existing full download (metadata.csv distinguishes it from an
        # apo_desolv-only download, which has no metadata.csv)
        if (target_dir / 'metadata.csv').is_file() and (
            target_dir / 'aligned_files'
        ).is_dir():
            mrich.print('Using existing hit data download', target_dir)
            self._hit_data_path = target_dir
            return target_dir

        path = DownloadService.download_target(
            target_name=target_name,
            target_access_string=project_name,
            proteins='',  # all observations (no poses exist yet to filter by)
            stack=stack,
            auth_token=auth_token,
            destination=destination,
            **{flag: True for flag in HIT_DATA_FLAGS},
        )

        self._hit_data_path = path
        return path

    def _ensure_apo_desolv_files(
        self, auth_token: str | None = None, stack: str = 'production'
    ) -> Path:
        """Ensure this target's apo-desolvated PDB files are available locally.

        Downloads the ``apo_desolv`` structures for this target's observations
        from Fragalysis (via :class:`.DownloadService`) and returns the path to
        the extracted directory. The download is performed at most once per
        instance and re-fetched each instance (overwriting any on-disk
        apo_desolv files) so a fresh instance always works with current data. If
        the full hit data was already downloaded this instance (via
        :meth:`._ensure_hit_data`), that is reused since it already includes
        fresh apo_desolv files.

        Everything needed for the request is taken from this animal: the target
        name and project (target access string) from :attr:`.target`, and the
        observation shortcodes from the ``pose_alias`` of this target's poses.

        .. note::
            This is a HIPPO-level helper, intended to be called from user-facing
            :class:`.HIPPO` methods the first time PDB files are needed. It's not
            expected to be called from the components or services layer. This
            method won't be necessary once HIPPO functions as a web service as
            intended.

        :param auth_token: optional Fragalysis ``sessionid``; otherwise the
            ``FRAGALYSIS_AUTH_TOKEN`` environment variable is used
        :param stack: Fragalysis stack to download from, a key into
            :data:`.STACK_URLS` (e.g. ``'production'``, ``'staging'``,
            ``'localhost'``); defaults to ``'production'``
        :returns: path to the extracted download directory
        """

        # already resolved during this session?
        if self._apo_desolv_path is not None and self._apo_desolv_path.exists():
            return self._apo_desolv_path

        # the full hit data downloaded this instance already includes fresh
        # apo_desolv files, so reuse it instead of re-downloading the subset
        if self._hit_data_path is not None and self._hit_data_path.exists():
            self._apo_desolv_path = self._hit_data_path
            return self._apo_desolv_path

        target_name = self._target.target_name
        project_name = self._target.project.project_name

        # downloads are laid out as data/downloads/<project_name>/<target_name>;
        # DownloadService extracts into destination/<target_name>, so we pass
        # data/downloads/<project_name> as the destination
        destination = DOWNLOADS_DIR / project_name

        # observation shortcodes to request, from the database
        proteins = list(
            PoseModel.objects.filter(target=self._target)
            .exclude(pose_alias__isnull=True)
            .exclude(pose_alias='')
            .values_list('pose_alias', flat=True)
            .distinct()
        )
        if not proteins:
            raise ValueError(
                f'No pose aliases found for target {target_name!r}; '
                'cannot determine which structures to download'
            )

        path = DownloadService.download_target(
            target_name=target_name,
            target_access_string=project_name,
            proteins=','.join(proteins),
            stack=stack,
            auth_token=auth_token,
            destination=destination,
        )

        self._apo_desolv_path = path
        self._apo_desolv_downloaded_at = datetime.now()
        return path

    def add_hits(
        self,
        *,
        metadata_csv: str | Path | None = None,
        aligned_directory: str | Path | None = None,
        auth_token: str | None = None,
        stack: str = 'production',
        tags: list | None = None,
        pose_methods: list[str] | None = None,
        skip: list | None = None,
        check_rmsd: bool = False,
        rmsd_threshold: float = 1.0,
        # debug: bool = False,
        # load_pose_mols: bool = False,
    ) -> pd.DataFrame:
        """Crystallographic hits from a Fragalysis download or XChemAlign alignment.

        Provide both `metadata_csv` and `aligned_directory` to load existing
        local data (for a Fragalysis download these point to the `metadata.csv`
        and `aligned_files` at the root of the extracted download; for an
        XChemAlign dataset `aligned_directory` points to the `aligned_files`).
        Omit both to download this target's data from the Fragalysis stack first
        (see :meth:`._ensure_hit_data`).

        :param metadata_csv: Path to the metadata.csv (omit to download)
        :param aligned_directory: Path to the aligned_files directory
            (omit to download)
        :param auth_token: optional Fragalysis ``sessionid`` for the download
            (otherwise ``FRAGALYSIS_AUTH_TOKEN`` is used)
        :param stack: Fragalysis stack to download from (default ``'production'``)
        :param skip: optional list of observation names to skip
        :returns: a DataFrame of metadata

        """

        ### Resolve the data source
        # Path-driven: provide both metadata_csv and aligned_directory to load
        # existing local data, or omit both to download the target's data from
        # the Fragalysis stack (always Fragalysis-type).
        if metadata_csv is None and aligned_directory is None:
            hit_dir = self._ensure_hit_data(auth_token=auth_token, stack=stack)
            aligned_directory = hit_dir / 'aligned_files'
            metadata_csv = hit_dir / 'metadata.csv'
        elif metadata_csv is None or aligned_directory is None:
            raise ValueError(
                'Provide both metadata_csv and aligned_directory to use existing '
                'data, or neither to download from the stack.'
            )

        skip = skip or []
        tags = tags or []
        pose_methods = pose_methods or DEFAULT_POSE_METHODS

        if not isinstance(aligned_directory, Path):
            aligned_directory = Path(aligned_directory)

        mrich.var('aligned_directory', aligned_directory)

        ### Validate inputs early with clear messages. A wrong/mismatched target
        # name usually yields an aligned_directory (often derived from the target
        # name) that doesn't exist or has no recognizable observation
        # subdirectories; without these checks that surfaces later as a confusing
        # "Unexpected mixed data format" assertion. We rely only on the aligned
        # data structure here -- not on the directory name, and not on the
        # presence of metadata (which is optional, e.g. for XChemAlign data).
        target_name = self.target.target_name
        if not aligned_directory.is_dir():
            raise NotADirectoryError(
                f'aligned_directory not found: {aligned_directory}. Check the path '
                f'matches the data for target {target_name!r}.'
            )

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

        subdirs = [p for p in aligned_directory.glob('*') if p.is_dir()]
        if not subdirs:
            raise ValueError(
                f'No observation subdirectories found in {aligned_directory}. Is the '
                'path correct and the download extracted? A wrong target name (here '
                f'{target_name!r}) often points add_hits at an empty/missing directory.'
            )

        SUBDIR_PATTERN_FRAGALYSIS = re.compile(r'^.*\d{4}[a-z]$')
        SUBDIR_PATTERN_XCA = re.compile(r'^.*-.\d{4}$')

        fragalysis_subdirs_present = any(
            SUBDIR_PATTERN_FRAGALYSIS.match(subdir.name) for subdir in subdirs
        )
        xca_subdirs_present = any(
            SUBDIR_PATTERN_XCA.match(subdir.name) for subdir in subdirs
        )

        # distinguish the two failure modes the old XOR assertion conflated
        if fragalysis_subdirs_present and xca_subdirs_present:
            raise ValueError(
                'Mixed Fragalysis and XChemAlign observation directories in '
                f'{aligned_directory}; expected a single consistent format.'
            )
        if not (fragalysis_subdirs_present or xca_subdirs_present):
            examples = ', '.join(p.name for p in subdirs[:3])
            raise ValueError(
                'Could not recognise any Fragalysis or XChemAlign observation '
                f'directories in {aligned_directory} (e.g. {examples}). Check that '
                f'the data matches target {target_name!r}.'
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

        pose_method_objs = []
        for name in pose_methods:
            obj = PoseMethodModel.objects.filter(pose_method_name=name).first()
            if obj is None:
                raise ValueError(
                    f"Pose method '{name}' not found. "
                    "Call register_pose_method() first."
                )
            pose_method_objs.append(obj)

        try:
            with transaction.atomic():
                result: IngestionBatchResult = IngestionService.ingest_filesystem(
                    root_path=aligned_directory,
                    target=self.target,
                    skip_records=skip,
                    compound_tag_list=tags,
                    metadata_file=metadata_csv,
                    pose_methods=pose_method_objs,
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
        enumeration_method: tuple[str, str] | None = None,
        pose_method: tuple[str, str] | None = None,
        score_cols: list[str] | None = None,
        scoring_methods: list[tuple[str, str]] | None = None,
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

        if name_col is None:
            raise ValueError(
                "name_col cannot be None. Provide the SDF column name that contains pose identifiers."
            )

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

        enumeration_method_obj = None
        if enumeration_method is not None:
            name, version = enumeration_method
            enumeration_method_obj = EnumerationMethodModel.objects.filter(
                enum_name=name, enum_version=version
            ).first()
            if enumeration_method_obj is None:
                raise ValueError(
                    f"Enumeration method '{name}' v{version} not found. "
                    "Call register_enumeration_method() first."
                )

        pose_method_obj = None
        if pose_method is not None:
            name, version = pose_method
            pose_method_obj = PoseMethodModel.objects.filter(
                pose_method_name=name, pose_method_version=version
            ).first()
            if pose_method_obj is None:
                raise ValueError(
                    f"Pose method '{name}' v{version} not found. "
                    "Call register_pose_method() first."
                )

        score_method_map = {}
        if score_cols and scoring_methods:
            if len(score_cols) != len(scoring_methods):
                raise ValueError('score_cols and scoring_methods must be the same length')
            for col, (method_name, method_version) in zip(score_cols, scoring_methods):
                obj = ScoringMethodModel.objects.filter(
                    method_name=method_name, method_version=method_version
                ).first()
                if obj is None:
                    raise ValueError(
                        f"Scoring method '{method_name}' v{method_version} not found. "
                        "Call register_scoring_method() first."
                    )
                score_method_map[col] = obj

        warn = make_warn_once_per_key()

        try:
            with transaction.atomic():
                result: IngestionBatchResult = IngestionService.ingest_sdf(
                    file_path=path,
                    target=self.target,
                    compound_tag_list=compound_tags,
                    pose_tag_list=pose_tags,
                    enumeration_method_obj=enumeration_method_obj,
                    pose_method_obj=pose_method_obj,
                    score_method_map=score_method_map,
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

    def add_enamine_real_routes(
        self,
        csv_path: str | Path,
        check_chemistry: bool = True,
        register_routes: bool = True,
    ) -> pd.DataFrame:
        """Add synthesis routes from an Enamine REAL CSV export"""

        try:
            with transaction.atomic():
                result = IngestionService.ingest_enamine_real_routes(
                    csv_path=csv_path,
                    do_check_chemistry=check_chemistry,
                    register_routes=register_routes,
                )
        except Exception as exc:
            logger.error(exc, exc_info=True)
            raise Exception from exc

        return result

    def prune_duplicate_routes(self) -> int:
        """Remove duplicate routes from the database"""
        return RouteService.prune_duplicate_routes()

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

    def set_derivative_subsites(self) -> None:
        """Propagate subsite assignments from inspiration poses to their derivatives."""
        SubsiteService.set_derivative_subsites()

    def register_enumeration_method(self, name: str, version: str, description: str = ''):
        """Register an enumeration method, or retrieve it if already registered."""
        return MethodService.register_enumeration_method(name, version, description)

    def register_pose_method(self, name: str, version: str, description: str = ''):
        """Register a pose method, or retrieve it if already registered."""
        return MethodService.register_pose_method(name, version, description)

    def register_scoring_method(self, name: str, version: str, description: str = ''):
        """Register a scoring method, or retrieve it if already registered."""
        return MethodService.register_scoring_method(name, version, description)

    @property
    def enumeration_methods(self):
        """All registered enumeration methods."""
        return MethodService.get_enumeration_methods()

    @property
    def pose_methods(self):
        """All registered pose methods."""
        return MethodService.get_pose_methods()

    @property
    def scoring_methods(self):
        """All registered scoring methods."""
        return MethodService.get_scoring_methods()
