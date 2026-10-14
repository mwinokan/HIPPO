"""Main animal class for HIPPO"""

import logging
import re
from datetime import datetime
from enum import Enum
from pathlib import Path

import mrich
import pandas as pd
from django.db import transaction

from .client import (
    GeneratorManager,
    IngredientManager,
    RecipeManager,
    RouteManager,
    ScorerManager,
)
from .models import (
    CompoundModel,
    PoseMethodModel,
    PoseModel,
    Project,
    RouteModel,
    TargetModel,
)
from .services.compound import CompoundTagService
from .services.download import DownloadService
from .services.ingestion import IngestionBatchResult, IngestionService
from .services.method import MethodService
from .services.pose import PoseTagService
from .services.quote import QuoteService
from .services.reaction import ReactionService
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

# When True, HIPPO.__init__ downloads this target's apo_desolv protein PDBs from
# Fragalysis (see HIPPO._ensure_apo_desolv_files). Toggle this module-level flag
# to enable/disable download-on-init -- deliberately NOT read from the
# environment. Currently False because the Fragalysis download/auth services are
# down for maintenance; set True to re-enable.
DOWNLOAD_APO_DESOLV_ON_INIT = False


class HIPPO:
    """Entry-point class of the xchem-hippo package.

    Update: this is atm not being called directly by the user.
    """

    def __init__(
        self,
        target_name: str,
        target_access_string: str,
        *,
        stack: str = 'production',
        auth_token: str | None = None,
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

        # Fragalysis stack / auth used for any downloads triggered by this
        # instance (see _ensure_hit_data / _ensure_apo_desolv_files).
        self._stack = stack
        self._auth_token = auth_token

        # Download state (see _ensure_hit_data / _ensure_apo_desolv_files). The
        # full hit data persists on disk and is reused across sessions; the
        # apo_desolv subset is re-downloaded once per instance to stay fresh.
        self._hit_data_path: Path | None = None
        self._apo_desolv_path: Path | None = None
        self._apo_desolv_downloaded_at: datetime | None = None

        # Optionally download this target's apo_desolv protein PDBs on init,
        # gated only by the DOWNLOAD_APO_DESOLV_ON_INIT flag (no longer requires a
        # previously downloaded aligned_files directory to be present).
        if DOWNLOAD_APO_DESOLV_ON_INIT:
            try:
                self._ensure_apo_desolv_files(
                    auth_token=self._auth_token, stack=self._stack
                )
            except ValueError as e:
                # no poses in the DB yet -> can't tell which structures to fetch
                mrich.warning(f'Skipping apo_desolv refresh on init: {e}')

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
    def reactants(self) -> CompoundSet:
        """Compounds that are reactants of a reaction and not a product of any
        (leaf reactants / purchasable building blocks)."""
        return CompoundSet(list(ReactionService.reactant_compound_ids()))

    @property
    def recipes(self) -> RecipeManager:
        """Client-side accessor for building recipes (see :class:`.RecipeManager`)."""
        return RecipeManager(self)

    @property
    def ingredients(self) -> IngredientManager:
        """Client-side accessor for building IngredientSets (see
        :class:`.IngredientManager`)."""
        return IngredientManager(self)

    @property
    def routes(self) -> RouteManager:
        """Client-side accessor for building RouteSets (see
        :class:`.RouteManager`)."""
        return RouteManager(self)

    @property
    def generators(self) -> GeneratorManager:
        """Client-side accessor for the random recipe/selection generators (see
        :class:`.GeneratorManager`)."""
        return GeneratorManager(self)

    @property
    def scorers(self) -> ScorerManager:
        """Client-side accessor for recipe scoring (see :class:`.ScorerManager`)."""
        return ScorerManager(self)

    @property
    def num_poses(self) -> int:
        """Total number of Poses in the Database"""
        return self.poses.count()

    def quote_compounds(
        self, compounds: 'CompoundSet | None' = None
    ) -> tuple[CompoundSet, CompoundSet]:
        """Report which compounds have catalogue quotes.

        :param compounds: optional :class:`.CompoundSet` to restrict to; defaults
            to all compounds in the database
        :returns: ``(quoted, unquoted)`` :class:`.CompoundSet` objects
        """
        if compounds is None:
            compounds = self.compounds
        elif not isinstance(compounds, CompoundSet):
            raise TypeError(
                f'compounds must be a CompoundSet or None, got {type(compounds)}'
            )

        quoted_ids, unquoted_ids = QuoteService.partition_quoted(compounds.ids)

        mrich.var('#quoted compounds', len(quoted_ids))
        mrich.var('#unquoted compounds', len(unquoted_ids))

        return CompoundSet(list(quoted_ids)), CompoundSet(list(unquoted_ids))

    def quote_reactants(self) -> tuple[CompoundSet, CompoundSet]:
        """Report which reactant compounds have catalogue quotes.

        Convenience wrapper around :meth:`.quote_compounds` restricted to the
        animal's reactants (see :attr:`.reactants`).

        :returns: ``(quoted, unquoted)`` :class:`.CompoundSet` objects
        """
        return self.quote_compounds(self.reactants)

    def plot_interaction_punchcard(
        self, poses: 'PoseSet | None' = None, *, subtitle=None, opacity=1.0, **kwargs
    ):
        """Plot an interaction punch-card for a :class:`.PoseSet` (default: all of
        this target's poses). See :func:`.plotting.plot_interaction_punchcard`."""
        from .plotting import plot_interaction_punchcard

        return plot_interaction_punchcard(
            self.poses if poses is None else poses,
            title_prefix=self._target.target_name,
            subtitle=subtitle,
            opacity=opacity,
            **kwargs,
        )

    def _ensure_hit_data(
        self, auth_token: str | None = None, stack: str = 'production'
    ) -> Path:
        """Ensure the target's full Fragalysis hit data is downloaded locally.

        Downloads all observations via :class:`.DownloadService` to
        ``data/downloads/<project>/<target>`` and returns that path. The data
        persists: an existing download (detected by ``metadata.csv``) is reused.

        :param auth_token: Fragalysis ``sessionid`` (else ``FRAGALYSIS_AUTH_TOKEN``)
        :param stack: Fragalysis stack to download from (default ``'production'``)
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
        """Ensure the target's apo-desolvated protein PDBs are downloaded locally.

        Downloads the ``apo_desolv`` structures for this target's poses (by
        ``pose_alias``) via :class:`.DownloadService`, once per instance. Reuses the
        full hit data if it was already downloaded this instance.

        :param auth_token: Fragalysis ``sessionid`` (else ``FRAGALYSIS_AUTH_TOKEN``)
        :param stack: Fragalysis stack to download from (default ``'production'``)
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
        stack: str | None = None,
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
        :param tags: compound tags to assign to every loaded hit
        :returns: a DataFrame of metadata
        :raises MissingTagError: if any tag is not registered. Ingestion never
            creates tags (see :meth:`register_compound_tag` /
            :meth:`register_pose_tag`), so this covers ``tags`` *and* every
            curated ``metadata.csv`` column that is used as a pose tag

        """

        # fall back to the stack/auth configured at instantiation
        if stack is None:
            stack = self._stack
        if auth_token is None:
            auth_token = self._auth_token

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
                    'Call register_pose_method() first.'
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
        max_workers: int | None = None,
        batch_size: int | None = None,
        chunk_size: int | None = None,
        single_transaction: bool = False,
        check_rmsd: bool = False,
        rmsd_threshold: float = 1.0,
    ) -> None:
        """Add posed virtual hits from an SDF into the database.

        Ingestion is set-based: the SDF is resolved to compounds, poses,
        inspirations, tags and scores in a fixed number of statements rather than a
        few dozen per record.

        All non-name columns are added to the pose metadata. N.B. separate .mol
        files are not created -- the molecule binary is stored in the database and
        fake paths are recorded.

        :param path: Path to the SDF
        :param reference: Optional single reference :class:`.PoseModel` used as the
            protein conformation for all poses, defaults to ``None``
        :param reference_col: Column containing reference :class:`.PoseModel`
            aliases or IDs, resolved per record
        :param inspirations: Optional :class:`.PoseSet` or list of IDs assigned as
            inspirations to every inserted pose, defaults to ``None``
        :param inspiration_col: Column containing per-record inspiration
            :class:`.PoseModel` names or IDs, defaults to ``"ref_mols"``
        :param inspiration_map: Optional mapping between inspiration strings found
            in ``inspiration_col`` and :class:`.PoseModel` ids
        :param compound_tags: String tags to assign to all created compounds,
            defaults to ``None``. Tags must already be registered (see
            :meth:`register_compound_tag`); ingestion never creates tags
        :param pose_tags: String tags to assign to all created poses, defaults to
            ``None``. Must already exist, as for ``compound_tags``
        :param enumeration_method: ``(name, version)`` of a registered enumeration
            method to associate with every compound
        :param pose_method: ``(name, version)`` of a registered pose method to
            associate with every pose
        :param score_cols: Columns holding score values
        :param scoring_methods: ``(name, version)`` pairs positionally aligned with
            ``score_cols``
        :param mol_col: Column containing the ``rdkit.ROMol`` ligands, defaults to
            ``"ROMol"``
        :param name_col: Column containing the ligand name/alias, defaults to
            ``"ID"``
        :param convert_floats: Try to convert all values to ``float``, defaults to
            ``True``
        :param skip_equal_dict: Skip rows where
            ``any(row[key] == value for key, value in skip_equal_dict.items())``,
            defaults to ``None``
        :param skip_not_equal_dict: Skip rows where
            ``any(row[key] != value for key, value in skip_not_equal_dict.items())``,
            defaults to ``None``
        :param max_workers: Worker processes for registration hashing, defaulting to
            ``min(8, cpu_count())``. Pass ``1`` to force serial hashing.
        :param batch_size: Rows per INSERT/UPDATE statement, or ``None`` for the
            automatic per-model maximum. Bounds *statement* size, not memory.
        :param chunk_size: SDF records read and processed per pass, defaulting to
            :data:`DEFAULT_CHUNK_SIZE`. This is what bounds *memory*.
        :param single_transaction: Wrap the whole file in one transaction rather
            than committing each chunk as it completes
        :param check_rmsd: Skip a pose whose RMSD to an existing pose of the same
            compound is below ``rmsd_threshold``
        :param rmsd_threshold: RMSD below which two poses are the same, in
            Angstrom, defaults to ``1.0``
        :raises MissingTagError: if any tag is not in the tag vocabulary

        .. note::
           Registration hashing dominates CPU cost. Compounds already registered
           under the same SMILES are resolved by an indexed lookup and never
           hashed; the rest are hashed across ``max_workers`` processes.

        .. note::
           ``chunk_size`` and ``batch_size`` control different things and are not
           interchangeable:

           * ``chunk_size`` is how many SDF records are held in memory at once. The
             file is streamed, so peak memory is roughly ``chunk_size`` x 100 KB
             regardless of file size -- a 300,000-record SDF loads in the same
             footprint as a 5,000-record one.
           * ``batch_size`` is how many rows go into a single INSERT. It exists
             because PostgreSQL caps a statement at 65535 bind parameters, which an
             unchunked ``bulk_create`` would breach at around 4,400 poses. It
             defaults to the largest safe value per model and is clamped down if an
             explicit value would breach the cap.

           They interact only in that a chunk cannot produce more rows than it
           holds: with ``chunk_size`` below ``batch_size``, each chunk is a single
           statement anyway.

        .. note::
           Re-ingesting a scored SDF updates existing scores rather than failing on
           the ``pk_score_values`` primary key.

        .. note::
           By default each chunk is committed as it completes, so a failure part way
           through leaves earlier chunks in the database. Ingestion is re-runnable --
           existing compounds and poses are found rather than duplicated -- so the
           fix is to run the same file again. Pass ``single_transaction=True`` for
           all-or-nothing semantics, at the cost of a transaction (and its locks and
           WAL) living for the entire load.
        """
        if not isinstance(path, Path):
            path = Path(path)

        if name_col is None:
            raise ValueError(
                'name_col cannot be None. Provide the SDF column name that '
                'contains pose identifiers.'
            )

        skip_equal_dict = skip_equal_dict or {}
        skip_not_equal_dict = skip_not_equal_dict or {}

        mrich.debug(f'{path=}')

        compound_tags = compound_tags or []
        pose_tags = pose_tags or []

        if isinstance(inspirations, PoseSet):
            inspiration_list = list(inspirations.ids)
        elif isinstance(inspirations, list):
            inspiration_list = inspirations
        else:
            inspiration_list = []

        if reference and isinstance(reference, PoseModel):
            reference_id = reference.id
        else:
            reference_id = None

        if inspiration_map is None:
            inspiration_map = {}

        # Method names are resolved to rows by the service layer (MethodService /
        # ScoreService), not here: that is DB traversal, and keeping ORM objects out
        # of the client-facing call is what the eventual client/backend split needs.

        warn = make_warn_once_per_key()

        # NB: deliberately does *not* disable trg_score_values_refresh_pivoted_mv.
        # That trigger is FOR EACH STATEMENT and scores are written one statement
        # per chunk, so disabling it would buy little while requiring schema helpers
        # that are not present on every deployment and risking leaving it off.
        #
        # Transaction scope lives in the service: the file is streamed in chunks and
        # each is committed as it completes (or all of them together, under
        # single_transaction), so there is no outer atomic() here.
        try:
            result: IngestionBatchResult = IngestionService.ingest_sdf(
                file_path=path,
                target=self.target,
                compound_tag_list=compound_tags,
                pose_tag_list=pose_tags,
                enumeration_method=enumeration_method,
                pose_method=pose_method,
                score_cols=score_cols,
                scoring_methods=scoring_methods,
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
                max_workers=max_workers,
                batch_size=batch_size,
                chunk_size=chunk_size,
                single_transaction=single_transaction,
                check_rmsd=check_rmsd,
                rmsd_threshold=rmsd_threshold,
            )
        except Exception:
            # re-raise the original rather than `raise Exception from exc`: a bare
            # Exception discards the type and traceback of the real failure
            logger.exception('load_sdf failed for %s', path)
            raise

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

        return result

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
            ["syndirella_product", "syndirella_placed"]. Tags must already be
            registered (see :meth:`register_pose_tag`); ingestion never creates
            tags, so these defaults must be present in the database
        :param product_tags: Add these tags to all inserted product compounds,
            defaults to ["syndirella_product"]. Must already exist, as for
            ``pose_tags``
        :returns: annotated DataFrame
        :raises MissingTagError: if any tag is not in the tag vocabulary
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
                    product_tag_list=product_tags,
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

    def register_enumeration_method(
        self, name: str, version: str, description: str = ''
    ):
        """Register an enumeration method, or retrieve it if already registered.

        Adding a new method needs DesignDB admin rights (from your database
        login); without them this raises :class:`.InsufficientPrivilegeError`.
        """
        return MethodService.register_enumeration_method(name, version, description)

    def register_pose_method(self, name: str, version: str, description: str = ''):
        """Register a pose method, or retrieve it if already registered.

        Adding a new method needs DesignDB admin rights (from your database
        login); without them this raises :class:`.InsufficientPrivilegeError`.
        """
        return MethodService.register_pose_method(name, version, description)

    def register_scoring_method(self, name: str, version: str, description: str = ''):
        """Register a scoring method, or retrieve it if already registered.

        Adding a new method needs DesignDB admin rights (from your database
        login); without them this raises :class:`.InsufficientPrivilegeError`.
        """
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

    def register_compound_tag(
        self, name: str, description: str | None = None, note: str | None = None
    ):
        """Register a compound tag, or retrieve it if already registered.

        Adding a new tag needs DesignDB admin rights (from your database login);
        without them this raises :class:`.InsufficientPrivilegeError`.
        """
        return CompoundTagService.register_tag(name, description, note)

    def register_pose_tag(
        self, name: str, description: str | None = None, note: str | None = None
    ):
        """Register a pose tag, or retrieve it if already registered.

        Adding a new tag needs DesignDB admin rights (from your database login);
        without them this raises :class:`.InsufficientPrivilegeError`.
        """
        return PoseTagService.register_tag(name, description, note)

    @property
    def compound_tags(self):
        """All registered compound tags."""
        return CompoundTagService.get_tags()

    @property
    def pose_tags(self):
        """All registered pose tags."""
        return PoseTagService.get_tags()
