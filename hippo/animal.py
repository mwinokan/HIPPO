import pandas as pd
import numpy as np

from .cset import CompoundTable, IngredientSet, CompoundSet
from .pset import PoseTable, PoseSet
from .tags import TagTable
from .rset import ReactionTable, ReactionSet
from .iset import InteractionTable
from .compound import Compound
from .reaction import Reaction
from .target import Target
from .pose import Pose
from .pycule import Quoter

from .db import Database
from pathlib import Path

from tqdm import tqdm

from .tools import inchikey_from_smiles, sanitise_smiles, SanitisationError

from mlog import setup_logger

logger = setup_logger("HIPPO")

from rdkit.Chem import Mol


class HIPPO:
    """The :class:`.HIPPO` `animal` class. Instantiating a :class:`.HIPPO` object will create or link a :class:`.HIPPO` :class:`.Database`.

    ::

            from hippo import HIPPO
            animal = HIPPO(project_name, db_path)

    .. attention::

            In addition to this API reference please see the tutorial pages :doc:`getting_started` and :doc:`insert_elaborations`.

    :param project_name: give this :class:`.HIPPO` a name
    :param db_path: path where the :class:`.Database` will be stored
    :returns: :class:`.HIPPO` object
    """

    def __init__(
        self,
        name: str,
        db_path: str | Path,
        update_legacy: bool = False,
    ):

        logger.header("Creating HIPPO animal")

        self._name = name
        # self._target_name = target

        logger.var("name", name, dict(color="arg"))
        # logger.var('target', target, dict(color='arg'))

        if not isinstance(db_path, Path):
            db_path = Path(db_path)

        logger.var("db_path", db_path, dict(color="file"))

        self._db_path = db_path
        self._db = Database(self.db_path, animal=self, update_legacy=update_legacy)

        self._compounds = CompoundTable(self.db)
        self._poses = PoseTable(self.db)
        self._tags = TagTable(self.db)
        self._reactions = ReactionTable(self.db)
        # self._interactions = InteractionTable(self.db)

        ### in memory subsets
        self._reactants = None
        self._products = None
        self._intermediates = None
        self._bases = None
        self._elabs = None

        logger.success(f"Initialised animal {self}")

    ### FACTORIES

    ### PROPERTIES

    @property
    def name(self) -> str:
        """Returns the project name

        :returns: project name
        """
        return self._name

    @property
    def db_path(self) -> str:
        """Returns the database path"""
        return self._db_path

    @property
    def db(self) -> Database:
        """Returns the Database object"""
        return self._db

    @property
    def compounds(self) -> CompoundTable:
        """Access compounds in the Database"""
        return self._compounds

    @property
    def poses(self) -> PoseTable:
        """Access Poses in the Database"""
        return self._poses

    @property
    def reactions(self) -> ReactionTable:
        """Access Reactions in the Database"""
        return self._reactions

    @property
    def tags(self) -> TagTable:
        """Access Tags in the Database"""
        return self._tags

    @property
    def interactions(self) -> InteractionTable:
        """Access Interactions in the Database"""
        # return self._interactions
        from .iset import InteractionTable

        return InteractionTable(self.db)

    @property
    def num_compounds(self) -> int:
        """Total number of Compounds in the Database"""
        return len(self.compounds)

    @property
    def num_poses(self) -> int:
        """Total number of Poses in the Database"""
        return len(self.poses)

    @property
    def num_reactions(self) -> int:
        """Total number of Reactions in the Database"""
        return len(self.reactions)

    @property
    def num_tags(self) -> int:
        """Number of unique Tags in the Database"""
        return len(self.tags.unique)

    @property
    def targets(self) -> list[Target]:
        """Access Targets in the Database"""
        target_ids = self.db.select(table="target", query="target_id", multiple=True)
        return [self.db.get_target(id=q) for q, in target_ids]

    @property
    def reactants(self) -> CompoundSet:
        """Returns all compounds that are reactants for at least one :class:`.Reaction` (and not products of others)"""
        if (
            self._reactants is None
            or self._reactants["total_changes"] != self.db.total_changes
        ):
            self._reactants = dict(
                set=self.compounds.reactants, total_changes=self.db.total_changes
            )
        return self._reactants["set"]

    @property
    def products(self) -> CompoundSet:
        """Returns all compounds that are products of at least one :class:`.Reaction` (and not reactants of others)"""
        if (
            self._products is None
            or self._products["total_changes"] != self.db.total_changes
        ):
            self._products = dict(
                set=self.compounds.products, total_changes=self.db.total_changes
            )
        return self._products["set"]

    @property
    def intermediates(self) -> CompoundSet:
        """Returns all compounds that are products and reactants of :class:`.Reaction`"""
        if (
            self._intermediates is None
            or self._intermediates["total_changes"] != self.db.total_changes
        ):
            self._intermediates = dict(
                set=self.compounds.intermediates, total_changes=self.db.total_changes
            )
        return self._intermediates["set"]

    @property
    def num_reactants(self) -> int:
        """Returns the number of reactants (see :meth:`reactants`)"""
        return len(self.reactants)

    @property
    def num_intermediates(self) -> int:
        """Returns the number of intermediates (see :meth:`intermediates`)"""
        return len(self.intermediates)

    @property
    def num_products(self) -> int:
        """Returns the number of products (see :meth:`products`)"""
        return len(self.products)

    @property
    def elabs(self) -> CompoundSet:
        """Returns compounds that are an based on another"""
        if self._elabs is None or self._elabs["total_changes"] != self.db.total_changes:
            self._elabs = dict(
                set=self.compounds.elabs, total_changes=self.db.total_changes
            )
        return self._elabs["set"]

    @property
    def bases(self) -> CompoundSet:
        """Returns compounds that are the basis for one or more elaborations"""
        if self._bases is None or self._bases["total_changes"] != self.db.total_changes:
            self._bases = dict(
                set=self.compounds.bases, total_changes=self.db.total_changes
            )
        return self._bases["set"]

    @property
    def num_elabs(self) -> int:
        """Number of compounds that are an elaboration of an existing base"""
        return len(self.elabs)

    @property
    def num_bases(self) -> int:
        """Number of compounds that are the basis for elaborations"""
        return len(self.bases)

    ### BULK INSERTION

    def add_hits(
        self,
        target_name: str,
        metadata_csv: str | Path,
        aligned_directory: str | Path,
        skip: list | None = None,
        debug: bool = False,
    ) -> pd.DataFrame:
        """Load in crystallographic hits downloaded from Fragalysis.

        :param target_name: Name of this protein :class:`.Target`
        :param metadata_csv: Path to the metadata.csv from the Fragalysis download
        :param aligned_directory: Path to the aligned_files directory from the Fragalysis download
        :param skip: optional list of observation names to skip
        :param debug: bool:  (Default value = False)
        :returns: a DataFrame of metadata

        """

        import molparse as mp
        from rdkit.Chem import PandasTools
        from .tools import remove_other_ligands

        if not isinstance(aligned_directory, Path):
            aligned_directory = Path(aligned_directory)

        # create the target
        target = self.register_target(name=target_name)

        skip = skip or []

        logger.var("aligned_directory", aligned_directory)

        count_directories_tried = 0
        count_compound_registered = 0
        count_poses_registered = 0

        meta_df = pd.read_csv(metadata_csv)
        generated_tag_cols = [
            "ConformerSites alias",
            "CanonSites alias",
            "CrystalformSites alias",
            "Quatassemblies alias",
            "Crystalforms alias",
        ]
        curated_tag_cols = [
            c
            for c in meta_df.columns
            if c
            not in ["Code", "Long code", "Compound code", "Smiles", "Downloaded"]
            + generated_tag_cols
        ]

        logger.var("curated_tag_cols", curated_tag_cols)

        for path in tqdm(list(aligned_directory.iterdir())):

            if not path.is_dir():
                continue

            if path.name in skip:
                continue

            count_directories_tried += 1

            sdfs = list(path.glob("*[0-9][0-9][0-9][0-9][a-z].sdf"))

            assert len(sdfs) == 1, (path, sdfs)

            pdbs = [
                p
                for p in path.glob("*.pdb")
                if "_ligand" not in p.name
                and "_apo" not in p.name
                and "_hippo" not in p.name
            ]

            assert len(pdbs) == 1, (path, pdbs)

            # load the SDF
            df = PandasTools.LoadSDF(
                str(sdfs[0]), molColName="ROMol", idName="ID", strictParsing=True
            )

            # extract fields
            observation_shortname = path.name.replace(".sdf", "")
            observation_longname = df.ID[0]
            mol = df.ROMol[0]
            lig_res_number = int(observation_longname.split("-")[1].split("_")[2])
            lig_chain = observation_longname.split("-")[1].split("_")[1]

            # smiles
            smiles = mp.rdkit.mol_to_smiles(mol)
            smiles = sanitise_smiles(smiles, verbosity=debug)

            # parse PDB
            sys = mp.parse(pdbs[0], verbosity=debug)

            if len(sys.residues["LIG"]) > 1 or any(
                r.contains_alternative_sites for r in sys.residues["LIG"]
            ):
                # create the single ligand bound pdb
                sys = remove_other_ligands(sys, lig_res_number, lig_chain)
                sys.prune_alternative_sites("A", verbosity=0)
                pose_path = str(pdbs[0].resolve()).replace(".pdb", "_hippo.pdb")
                mp.write(pose_path, sys, shift_name=True, verbosity=debug)
            else:
                pose_path = str(pdbs[0].resolve())

            # create the molecule / pose
            compound_id = self.db.insert_compound(
                # name=crystal_name,
                smiles=smiles,
                tags=["hits"],
                warn_duplicate=debug,
            )

            if not compound_id:

                inchikey = inchikey_from_smiles(smiles)
                compound = self.compounds[inchikey]

                if not compound:
                    logger.error(
                        "Compound exists in database but could not be found by inchikey"
                    )
                    logger.var("smiles", smiles)
                    logger.var("inchikey", inchikey)
                    logger.var("observation_shortname", observation_shortname)
                    raise Exception

            else:
                count_compound_registered += 1
                compound = self.compounds[compound_id]

            # metadata

            meta_row = meta_df[meta_df["Code"] == observation_shortname]
            if not len(meta_row):
                assert observation_longname
                meta_row = meta_df[meta_df["Long code"] == observation_longname]

            # logger.var('observation_longname',observation_longname)

            assert len(meta_row)

            tags = ["hits"]
            for tag in curated_tag_cols:
                if meta_row[tag].values[0]:
                    tags.append(tag)

            metadata = {"observation_longname": observation_longname}
            for tag in generated_tag_cols:
                metadata[tag] = meta_row[tag].values[0]

            pose_id = self.db.insert_pose(
                compound=compound,
                alias=observation_shortname,
                target=target.id,
                path=pose_path,
                tags=tags,
                metadata=metadata,
            )

            if pose_id:
                count_poses_registered += 1

        logger.var("#directories parsed", count_directories_tried)
        logger.var("#compounds registered", count_compound_registered)
        logger.var("#poses registered", count_poses_registered)

    def add_compounds(
        self,
        target: str,
        sdf_path: str | Path,
        *,
        reference: int | Pose | None = None,
        tags: None | list[str] = None,
        output_directory: str | Path | None = None,
        mol_col: str = "ROMol",
        name_col: str = "ID",
        inspiration_col: str = "ref_mols",
        inspiration_map: None | dict = None,
        skip_first: bool = False,
        convert_floats: bool = True,
        stop_after: int | None = None,
        skip_equal_dict: dict | None = None,
        skip_not_equal_dict: dict | None = None,
        check_pose_RMSD: bool = False,
        pose_RMSD_tolerance: float = 1.0,
    ) -> None:
        """Add posed virtual hits from an SDF into the database.

        :param target: Name of the protein :class:`.Target`
        :param sdf_path: Path to the SDF
        :param reference: Reference :class:`.Pose` to use as the protein conformation for all poses, defaults to ``None``
        :param tags: Tags to assign to all created compounds and poses, defaults to ``None``
        :param output_directory: Specify the path where individual ligand .mol files will be created, defaults to ``None`` where the name of the SDF is used.
        :param mol_col: Name of the column containing the ``rdkit.ROMol`` ligands, defaults to ``"ROMol"``
        :param name_col: Name of the column containing the ligand name/alias, defaults to ``"ID"``
        :param inspiration_col: Name of the column containing the list of inspiration pose names, defaults to ``"ref_mols"``
        :param inspiration_map: Dictionary or callable mapping between inspiration strings found in ``inspiration_col`` and :class:`.Pose` ids, defaults to None
        :param skip_first: Skip first row, defaults to ``False``
        :param convert_floats: Try to convert all values to ``float``, defaults to ``True``
        :param stop_after: Stop after given number of rows, defaults to ``None``
        :param skip_equal_dict: Skip rows where ``any(row[key] == value for key, value in skip_equal_dict.items())``, defaults to ``None``
        :param skip_not_equal_dict: Skip rows where ``any(row[key] != value for key, value in skip_not_equal_dict.items())``, defaults to ``None``
        :param check_pose_RMSD: Check :class:`.Pose` to previously registered poses and skip if below ``pose_RMSD_tolerance``, defaults to ``False``
        :param pose_RMSD_tolerance: Tolerance for ``check_pose_RMSD`` in Angstrom, defaults to ``1.0``

        All non-name columns added to the Pose metadata.
        """

        if not isinstance(sdf_path, Path):
            sdf_path = Path(sdf_path)

        logger.debug(f"{sdf_path=}")

        tags = tags or []

        from rdkit.Chem import PandasTools, MolToMolFile, MolFromMolFile
        from molparse.rdkit import mol_to_smiles, mol_to_pdb_block
        from numpy import isnan
        from pandas import read_pickle

        if sdf_path.name.endswith(".sdf"):
            df = PandasTools.LoadSDF(str(sdf_path.resolve()))
        else:
            df = read_pickle(sdf_path)

        df_columns = list(df.columns)

        target = self.register_target(target).name

        assert mol_col in df_columns, f"{mol_col=} not in {df_columns}"
        assert name_col in df_columns, f"{name_col=} not in {df_columns}"
        assert inspiration_col in df_columns, f"{inspiration_col=} not in {df_columns}"

        df_columns.pop(df_columns.index(mol_col))
        df_columns.pop(df_columns.index(name_col))
        df_columns.pop(df_columns.index(inspiration_col))

        if not output_directory:
            import os

            output_directory = str(sdf_path.name).removesuffix(".sdf")

        output_directory = Path(output_directory)
        if not output_directory.exists:
            logger.writing(f"Creating output directory {output_directory}")
            os.system(f"mkdir -p {output_directory}")

        n_poses = self.num_poses
        n_comps = self.num_compounds

        for i, row in tqdm(df.iterrows()):

            name = row[name_col].strip()

            if name == "ver_1.2":
                logger.warning("Skipping Fragalysis header molecule")
                continue

            if skip_equal_dict and any(row[k] == v for k, v in skip_equal_dict.items()):
                continue

            if skip_not_equal_dict and any(
                row[k] != v for k, v in skip_not_equal_dict.items()
            ):
                continue

            mol = row[mol_col]

            if not name:
                name = f"pose_{i}"

            mol_path = output_directory / f"{name}.mol"

            MolToMolFile(mol, str(mol_path.resolve()))

            smiles = mol_to_smiles(mol)

            comp = self.register_compound(smiles=smiles, tags=tags)

            if not comp:
                logger.error(f"Could not register compound {i=}")
                continue

            # inspirations

            inspirations = []

            insp_str = row[inspiration_col]

            if isinstance(insp_str, str):
                insp_str = insp_str.removeprefix("[")
                insp_str = insp_str.removesuffix("]")
                insp_str = insp_str.replace("'", "")
                generator = insp_str.split(",")

            else:
                generator = insp_str

            for insp in generator:
                insp = insp.strip()

                if inspiration_map:
                    if isinstance(inspiration_map, dict) and insp in inspiration_map:
                        pose = inspiration_map[insp]
                        if pose:
                            inspirations.append(pose)
                    elif hasattr(inspiration_map, "__call__"):
                        pose = inspiration_map(insp)
                        if pose:
                            inspirations.append(pose)
                    else:
                        pose = self.poses[insp]
                        if pose:
                            inspirations.append(pose)
                        else:
                            logger.error(f"Could not find inspiration pose {insp}")
                            continue

                else:
                    pose = self.poses[insp]
                    if pose:
                        inspirations.append(pose)
                    else:
                        logger.error(f"Could not find inspiration pose {insp}")
                        continue

            # metadata
            metadata = {}

            for col in df_columns:
                value = row[col]

                if not isinstance(col, str):
                    if i == 0:
                        logger.warning(f"Skipping metadata from column={col}.")
                    continue

                if isinstance(value, float) and isnan(value):
                    continue

                if convert_floats:
                    try:
                        value = float(value)
                    except TypeError:
                        pass
                    except ValueError:
                        pass

                if not (isinstance(value, str) or isinstance(value, float)):
                    if i == 0:
                        logger.warning(f"Skipping metadata from column={col}.")
                    continue

                metadata[col] = value

            # print(metadata)

            pose = self.register_pose(
                alias=name,
                compound=comp,
                target=target,
                path=mol_path,
                metadata=metadata,
                inspirations=inspirations,
                tags=tags,
                reference=reference,
                check_RMSD=check_pose_RMSD,
                RMSD_tolerance=pose_RMSD_tolerance,
            )

            if stop_after and i + 1 >= stop_after:
                break

        if n := self.num_compounds - n_comps:
            f = logger.success
        else:
            f = logger.error

        f(f"Loaded {n} compounds from {sdf_path}")

        if n := self.num_poses - n_poses:
            f = logger.success
        else:
            f = logger.error

        f(f"Loaded {n} poses from {sdf_path}")

    def add_syndirella_bases(
        self,
        df_path: str | Path,
        *,
        tags: None | list[str] = None,
        reaction_yield_map: dict | None = None,
        check_chemistry: bool = True,
        return_reactions: bool = False,
        debug: bool = False,
    ) -> CompoundSet | tuple[CompoundSet, ReactionSet]:
        """
        Load Syndirella base compounds and poses from a pickled DataFrame

        :param df_path: Path to the pickled DataFrame or SDF.
        :param tags: list of tags to assign to compounds and poses, defaults to ``None``
        :param reaction_yield_map: dictionary mapping reaction type strings to their yield ratio, defaults to ``None``
        :param check_chemistry: check the reaction chemistry, defaults to ``True``
        :param return_reactions: Return :class:`.ReactionSet` of reactions, defaults to ``False``
        :param debug: Increase verbosity of output, defaults to ``False``
        :returns: :class:`.CompoundSet` of compounds if ``return_reactions = False``, else returns a tuple: ``(compounds, reactions)``
        """

        raise NotImplementedError

        from .chem import check_reaction_types, InvalidChemistryError

        tags = tags or []

        if isinstance(df_path, str):
            df_path = Path(df_path)

        from pickle import UnpicklingError

        if df_path.name.endswith(".csv"):
            raise NotImplementedError
        elif df_path.name.endswith(".sdf"):
            from rdkit.Chem import PandasTools

            df = PandasTools.LoadSDF(df_path)
        else:
            try:
                df = pd.read_pickle(df_path)
            except UnpicklingError:
                logger.error(
                    "Could not read DataFrame. Is it a valid sdf or pickled dataframe?"
                )
                return None

        if debug:
            logger.var("len(df)", len(df))

        reaction_cols = [c for c in df.columns if "_reaction" in c]

        product_ids = set()

        reactions = []

        # loop over rows

        for i, row in tqdm(df.iterrows()):

            # determine number of steps

            reaction_types = [
                r for c in reaction_cols if isinstance((r := row[c]), str)
            ]
            check_reaction_types(reaction_types)

            n_steps = len(reaction_types)

            if debug:
                logger.header(f"row={i}, {n_steps=}")

            try:

                # loop over each reaction step
                for j in range(n_steps):

                    j += 1

                    reactants = []

                    reactant_previous_product = row[f"{j}_r_previous_product"]

                    # reactant 1
                    if reactant_previous_product == 1:
                        reactant1_id = product_id
                        reactants.append(reactant1_id)
                    elif smiles := row[f"{j}_r1_smiles"]:

                        if not isinstance(smiles, str):
                            raise InvalidRowError(f"non-string {j}_r1_smiles")

                        reactant1_id, duplicate = self.register_compound(
                            smiles=smiles,
                            commit=False,
                            return_compound=False,
                            tags=tags,
                            base=None,
                            register_base_if_duplicate=False,
                            return_duplicate=True,
                        )

                        reactants.append(reactant1_id)

                    # reactant 2
                    if reactant_previous_product == 2:
                        reactant2_id = product_id
                        reactants.append(reactant2_id)
                    elif smiles := row[f"{j}_r2_smiles"]:

                        if not isinstance(smiles, str):
                            raise InvalidRowError(f"non-string {j}_r2_smiles")

                        reactant2_id, duplicate = self.register_compound(
                            smiles=smiles,
                            commit=False,
                            return_compound=False,
                            tags=tags,
                            base=None,
                            register_base_if_duplicate=False,
                            return_duplicate=True,
                        )

                        reactants.append(reactant2_id)

                    # product
                    if smiles := row[f"{j}_product_smiles"]:
                        if not isinstance(smiles, str):
                            raise InvalidRowError(f"non-string {j}_product_smiles")

                        if j != n_steps:
                            this_tags = tags

                        else:
                            this_tags = ["Syndirella base"] + tags

                        product_id, duplicate = self.register_compound(
                            smiles=smiles,
                            tags=this_tags,
                            commit=False,
                            return_compound=False,
                            base=None,
                            register_base_if_duplicate=False,
                            return_duplicate=True,
                        )

                        if j == n_steps:
                            product_ids.add(product_id)

                    # register the reaction
                    if reaction_yield_map:
                        product_yield = reaction_yield_map[row[f"{j}_reaction"]]
                    else:
                        product_yield = 1.0

                    try:

                        reaction = self.register_reaction(
                            reactants=reactants,
                            product=product_id,
                            type=row[f"{j}_reaction"],
                            commit=False,
                            product_yield=product_yield,
                            check_chemistry=check_chemistry,
                        )

                        reactions.append(reaction.id)

                        if debug:
                            logger.var(str(reaction), reaction.reaction_str)

                    except InvalidChemistryError as e:
                        skipped_reactions += 1

            except InvalidRowError as e:
                skipped_invalid_smiles += 1
                continue

        if return_reactions:
            return self.compounds[product_ids], self.reactions[reactions]
        else:
            return self.compounds[product_ids]

    def add_syndirella_elabs(
        self,
        df_path: str | Path,
        *,
        inspiration_map: dict | None = None,
        base_only: bool = False,
        tags: None | list[str] = None,
        reaction_yield_map: dict | None = None,
        require_truthy_bases: None | list[str] = None,
        # warn_nontruthy_bases: None | list[str] = None,
        require_truthy_elabs: None | list[str] = None,
        warn_nonzero_truthy_bases: None | list[str] = None,
        stop_after: None | int = None,
        check_chemistry: bool = True,
        base_designator: str = "scaffold",
        target: int | str = 1,
        allow_missing_inspirations: bool = False,
    ) -> None:
        """
        Load Syndirella elaboration compounds and poses from a pickled DataFrame

        :param df_path: Path to the pickled DataFrame
        :param inspiration_map: Dictionary or callable mapping between inspiration strings found in ``inspiration_col`` and :class:`.Pose` ids, defaults to None
        :param base_only: Only load the base compound, defaults to ``False``
        :param tags: list of tags to assign to compounds and poses, defaults to ``None``
        :param reaction_yield_map: dictionary mapping reaction type strings to their yield ratio, defaults to ``None``
        :param require_truthy_bases: List of columns that should be truthy to load the given base, defaults to ``['path_to_mol', 'intra_geometry_pass']``
        :param require_truthy_elabs: List of columns that should be truthy to load the given elaboration, defaults to ``['path_to_mol', 'intra_geometry_pass']``
        :param warn_nonzero_truthy_bases: List of columns that should have a truthy value for at least one of the base molecules, defaults to ``['path_to_mol', 'intra_geometry_pass']``
        :param stop_after: Stop after given number of rows, defaults to ``None``
        :param check_chemistry: check the reaction chemistry, defaults to ``True``
        """

        from .chem import check_reaction_types, InvalidChemistryError

        if inspiration_map is None:
            inspiration_map = self.db.create_metadata_id_map(
                table="pose", key="observation_longname"
            )

        tags = tags or []
        require_truthy_bases = require_truthy_bases or [
            "path_to_mol",
            "intra_geometry_pass",
        ]
        require_truthy_bases = require_truthy_bases or [
            "path_to_mol",
            "intra_geometry_pass",
        ]
        require_truthy_elabs = require_truthy_elabs or [
            "path_to_mol",
            "intra_geometry_pass",
        ]
        warn_nonzero_truthy_bases = warn_nonzero_truthy_bases or [
            "path_to_mol",
            "intra_geometry_pass",
        ]
        assert all([k in require_truthy_bases for k in warn_nonzero_truthy_bases])

        if isinstance(df_path, str):
            df_path = Path(df_path)
        df = pd.read_pickle(df_path)

        # work out number of reaction steps
        n_steps = max(
            [int(s.split("_")[0]) for s in df.columns if "_product_smiles" in s]
        )

        # check chemistries
        if check_chemistry:
            chemistries = set(
                sum(
                    [
                        df[df[f"{j+1}_reaction"].notnull()][f"{j+1}_reaction"].tolist()
                        for j in range(n_steps)
                    ],
                    [],
                )
            )
            logger.var("Present reactions", str(chemistries))
            check_reaction_types(chemistries)

        if f"{n_steps}_num_atom_diff" not in df.columns:
            logger.error(df_path)
            logger.error(f"{n_steps}_num_atom_diff not in columns:")
            print(df.columns)
            raise NotImplementedError

        base_df = df[df[f"{n_steps}_product_name"].str.contains(base_designator)]

        if not len(base_df):
            logger.error(f"No base/scaffold rows found in {df_path.name}")
            return df

        for key in warn_nonzero_truthy_bases:
            if not any(base_df[key].values):
                logger.warning(f"No bases have truthy {key}. Inserting them anyway")
                require_truthy_bases.pop(require_truthy_bases.index(key))

        base_id = None
        base_reactants = {}

        if base_only:
            generator = df.iterrows()
        else:
            generator = tqdm(df.iterrows())

        n_comps = len(self.compounds)
        n_poses = len(self.poses)
        n_reactions = len(self.reactions)

        skipped_smaller = 0
        skipped_reactions = 0
        skipped_invalid_smiles = 0

        try:

            for i, row in generator:

                path_to_mol = row.path_to_mol

                product_name = row[f"{n_steps}_product_name"]

                this_row_is_a_base = base_designator in product_name

                # skip entries that have non-truthy columns
                if this_row_is_a_base:

                    if any(not row[key] for key in require_truthy_bases):
                        logger.warning(
                            f"Skipping (row {i=}) which has {key}={row[key]}"
                        )
                        continue

                elif any(not row[key] for key in require_truthy_elabs):
                    continue

                if row[f"{n_steps}_num_atom_diff"] <= 0 and not this_row_is_a_base:
                    skipped_smaller += 1
                    continue

                if base_only and not this_row_is_a_base:
                    continue

                if base_only and base_id and not this_row_is_a_base:
                    break

                # logger.debug(f'{i=} {this_row_is_a_base=}')

                if this_row_is_a_base and skipped_smaller and not base_id:
                    logger.warning(
                        f"Skipped {skipped_smaller} elaborations that are smaller than the base compound"
                    )

                elabs_registered = set()

                try:

                    # loop over each reaction step
                    for j in range(n_steps):

                        j += 1

                        reactants = []

                        reactant_previous_product = row[f"{j}_r_previous_product"]

                        # reactant 1
                        if reactant_previous_product == 1:
                            reactant1_id = product_id
                            reactants.append(reactant1_id)
                        elif smiles := row[f"{j}_r1_smiles"]:

                            if not isinstance(smiles, str):
                                raise InvalidRowError(f"non-string {j}_r1_smiles")

                            try:
                                base = (
                                    base_reactants[j][1]
                                    if not this_row_is_a_base
                                    else None
                                )
                            except KeyError:
                                print(row)
                                logger.error(
                                    f"Expected base_reactants to contain data when {this_row_is_a_base=}"
                                )
                                raise

                            reactant1_id, duplicate = self.register_compound(
                                smiles=smiles,
                                commit=False,
                                return_compound=False,
                                tags=tags,
                                bases=[base],
                                register_base_if_duplicate=True,
                                return_duplicate=True,
                            )

                            reactants.append(reactant1_id)

                        # reactant 2
                        if reactant_previous_product == 2:
                            reactant2_id = product_id
                            reactants.append(reactant2_id)

                        elif smiles := row[f"{j}_r2_smiles"]:

                            if not isinstance(smiles, str):
                                raise InvalidRowError(f"non-string {j}_r2_smiles")

                            if this_row_is_a_base:
                                base = None
                            else:
                                base = base_reactants[j][2]
                                assert base

                            reactant2_id, duplicate = self.register_compound(
                                smiles=smiles,
                                commit=False,
                                return_compound=False,
                                tags=tags,
                                bases=[base],
                                register_base_if_duplicate=True,
                                return_duplicate=True,
                            )

                            reactants.append(reactant2_id)

                        else:
                            reactant2_id = None

                        # product
                        if smiles := row[f"{j}_product_smiles"]:
                            if not isinstance(smiles, str):
                                raise InvalidRowError(f"non-string {j}_product_smiles")

                            if j != n_steps:
                                this_tags = tags
                                base = None

                            elif this_row_is_a_base:
                                this_tags = ["Syndirella base"] + tags
                                base = None

                            else:
                                this_tags = tags
                                base = base_id

                            product_id, duplicate = self.register_compound(
                                smiles=smiles,
                                tags=this_tags,
                                commit=False,
                                return_compound=False,
                                bases=[base],
                                register_base_if_duplicate=True,
                                return_duplicate=True,
                            )

                            if not base_id:
                                base_reactants[j] = {
                                    1: reactant1_id,
                                    2: reactant2_id,
                                    "product": product_id,
                                }

                            if not base_id and j == n_steps and this_row_is_a_base:
                                base_id = product_id

                        # register the reaction
                        if reaction_yield_map:
                            product_yield = reaction_yield_map[row[f"{j}_reaction"]]
                        else:
                            product_yield = 1.0

                        try:
                            # print(reactants, product_id, row[f"{j}_reaction"])
                            self.register_reaction(
                                reactants=reactants,
                                product=product_id,
                                type=row[f"{j}_reaction"],
                                commit=False,
                                product_yield=product_yield,
                                check_chemistry=check_chemistry,
                            )
                        except InvalidChemistryError as e:
                            skipped_reactions += 1

                except InvalidRowError as e:
                    logger.error(f"Skipping invalid row {i=}: {e}")
                    skipped_invalid_smiles += 1
                    continue

                # pose metadata
                metadata = {}

                # inspirations
                inspirations = []
                for inspiration in row.regarded:
                    if inspiration_map:
                        match inspiration_map:
                            case dict():
                                try:
                                    inspiration = inspiration_map[inspiration]
                                except KeyError:
                                    logger.error(
                                        f"{inspiration=} not found in inspiration_map"
                                    )
                                    inspiration = None
                                    if not allow_missing_inspirations:
                                        from json import dumps

                                        print(dumps(inspiration_map, indent=2))
                                        raise
                            case _:
                                inspiration = inspiration_map(inspiration)

                    else:
                        # this is really expensive
                        inspiration = self.poses[inspiration]
                    if inspiration:
                        inspirations.append(inspiration)

                if path_to_mol:

                    path_to_mol = path_to_mol.replace("//", "/")

                    # register the pose
                    self.register_pose(
                        compound=product_id,
                        target=target,
                        path=path_to_mol,
                        inspirations=inspirations,
                        metadata=metadata,
                        tags=this_tags,
                        commit=False,
                        return_pose=False,
                        overwrite_metadata=True,
                        warn_duplicate=False,
                        energy_score=row["∆∆G"],
                        distance_score=row["comRMSD"],
                    )

                if stop_after and i > stop_after:
                    break

        except KeyboardInterrupt:
            logger.error("KeyboardInterrupt")

        self.db.commit()

        n_comps = len(self.compounds) - n_comps
        n_poses = len(self.poses) - n_poses
        n_reactions = len(self.reactions) - n_reactions

        if skipped_reactions:
            logger.warning(f"Skipped {skipped_reactions} invalid reactions")

        if skipped_invalid_smiles:
            logger.warning(f"Skipped {skipped_invalid_smiles} rows with NaN smiles")

        if n_comps:
            if not base_only:
                logger.success(f"Loaded {n_comps} new compounds from {df_path.name}")
        else:
            logger.warning(f"Loaded {n_comps} new compounds from {df_path.name}")

        if n_poses:
            if not base_only:
                logger.success(f"Loaded {n_poses} new poses from {df_path.name}")
        else:
            logger.warning(f"Loaded {n_poses} new poses from {df_path.name}")

        if n_reactions:
            if not base_only:
                logger.success(
                    f"Loaded {n_reactions} new reactions from {df_path.name}"
                )
        else:
            logger.warning(f"Loaded {n_reactions} new reactions from {df_path.name}")

        return df

    def add_enamine_quote(
        self,
        path: str | Path,
        *,
        orig_name_col: str = "Customer Code",
        # orig_name_col: str = 'Diamond ID (Molecule Name)',
        entry_col="Catalog ID",
        stop_after=None,
        orig_name_is_hippo_id: bool = False,
    ):
        """
        Load an Enamine quote provided as an excel file

        :param path: Path to the excel file
        :param orig_name_col: Column name of the original alias, defaults to 'Customer Code'
        :param entry_col: Column name of the catalogue ID/entry, defaults to 'Catalog ID'
        :param stop_after: Stop after given number of rows, defaults to ``None``
        :param orig_name_is_hippo_id: Set to ``True`` if ``orig_name_col`` is the original HIPPO :class:``hippo.compound.Compound`` ID, defaults to ``False``
        :returns: An :class:`.IngredientSet` of the quoted molecules
        """

        df = pd.read_excel(path)

        smiles_col = "SMILES"
        purity_col = "Purity, %"
        amount_col = "Amount, mg"
        catalogue_col = "Collection"
        lead_time_col = "Lead time"

        assert smiles_col in df.columns, "Unexpected Excel format (smiles_col)"
        assert orig_name_col in df.columns, "Unexpected Excel format (orig_name_col)"
        assert entry_col in df.columns, "Unexpected Excel format (entry_col)"
        assert purity_col in df.columns, "Unexpected Excel format (purity_col)"
        assert amount_col in df.columns, "Unexpected Excel format (amount_col)"
        assert catalogue_col in df.columns, "Unexpected Excel format (catalogue_col)"
        assert (
            "Price, EUR" in df.columns or "Price, USD" in df.columns
        ), "Unexpected Excel format (price_col)"
        assert lead_time_col in df.columns, "Unexpected Excel format (lead_time_col)"

        price_cols = [c for c in df.columns if c.startswith("Price")]
        assert len(price_cols) == 1
        price_col = price_cols[0]
        currency = price_col.split(", ")[-1]

        ingredients = IngredientSet(self.db)

        if len(df) > 100:
            generator = tqdm(df.iterrows(), total=len(df))
        else:
            generator = df.iterrows()

        for i, row in generator:
            smiles = row[smiles_col]

            if not isinstance(smiles, str):
                break

            compound = self.register_compound(smiles=smiles)

            if orig_name_is_hippo_id:
                try:
                    expected_id = int(row[orig_name_col])

                    if expected_id != compound.id:
                        logger.error("Compound registration mismatch:")
                        logger.var("expected_id", expected_id)
                        logger.var("new_id", compound.id)
                        logger.var(
                            "original_smiles", self.compounds[expected_id].smiles
                        )
                        logger.var("new_smiles", smiles)

                except ValueError:
                    pass

            if (catalogue := row[catalogue_col]) == "No starting material":
                continue

            if (price := row[price_col]) == 0.0:
                continue

            if not isinstance(row[lead_time_col], str):
                continue

            if "week" in row[lead_time_col]:
                lead_time = int(row[lead_time_col].split()[0].split("-")[-1]) * 5
            else:
                raise NotImplementedError

            quote_data = dict(
                compound=compound,
                supplier="Enamine",
                catalogue=catalogue,
                entry=row[entry_col],
                amount=row[amount_col],
                purity=row[purity_col] / 100,
                lead_time=lead_time,
                price=price,
                currency=currency,
                smiles=smiles,
            )

            q_id = self.db.insert_quote(**quote_data)

            ingredients.add(
                compound_id=compound.id,
                amount=row[amount_col],
                quoted_amount=row[amount_col],
                quote_id=q_id,
                supplier="Enamine",
                max_lead_time=None,
            )

            if stop_after and stop_after == i:
                break

        return ingredients

    def add_mcule_quote(
        self,
        path: str | Path,
    ):
        """
        Load an MCule quote provided as an excel file

        :param path: Path to the excel file
        :returns: An :class:`.IngredientSet` of the quoted molecules
        """

        ### get lead time from suppliers sheet

        sheet_name: str = "List of suppliers"
        df = pd.read_excel(path, sheet_name=sheet_name)

        supplier_col = "Supplier"
        lead_time_col = "Delivery time (working days)"

        assert supplier_col in df.columns, "Unexpected Excel format (supplier_col)"
        assert lead_time_col in df.columns, "Unexpected Excel format (lead_time_col)"

        lead_time_lookup = {
            row[supplier_col]: row[lead_time_col] for i, row in df.iterrows()
        }

        ### parse individual compound quotes

        sheet_name: str = "List of products"
        df = pd.read_excel(path, sheet_name=sheet_name)

        # return df

        smiles_col = "Quoted product SMILES"
        entry_col = "Query Mcule ID"
        purity_col = "Guaranteed purity (%)"
        amount_col = "Quoted Amount (mg)"
        catalogue_col = "Supplier"
        lead_time_col = "Lead time"
        price_col = "Product price (USD)"
        currency = "USD"

        assert smiles_col in df.columns, "Unexpected Excel format (smiles_col)"
        assert entry_col in df.columns, "Unexpected Excel format (entry_col)"
        assert purity_col in df.columns, "Unexpected Excel format (purity_col)"
        assert amount_col in df.columns, "Unexpected Excel format (amount_col)"
        assert catalogue_col in df.columns, "Unexpected Excel format (catalogue_col)"
        assert price_col in df.columns, "Unexpected Excel format (price_col)"

        ingredients = IngredientSet(self.db)

        for i, row in tqdm(df.iterrows()):
            smiles = row[smiles_col]

            if not isinstance(smiles, str):
                break

            compound = self.register_compound(smiles=smiles)

            # if (catalogue := row[catalogue_col]) == 'No starting material':
            #   continue

            catalogue = row[catalogue_col]
            lead_time = lead_time_lookup[catalogue]

            # if (price := row[price_col]) == 0.0:
            #   continue

            # if not isinstance(row[lead_time_col], str):
            # continue

            # if 'week' in row[lead_time_col]:
            #   lead_time = int(row[lead_time_col].split()[0].split('-')[-1])*5
            # else:
            #   raise NotImplementedError

            quote_data = dict(
                compound=compound,
                supplier="MCule",
                catalogue=catalogue,
                entry=row[entry_col],
                amount=row[amount_col],
                purity=row[purity_col] / 100,
                lead_time=lead_time,
                price=row[price_col],
                currency=currency,
                smiles=smiles,
            )

            q_id = self.db.insert_quote(**quote_data, commit=False)

            ingredients.add(
                compound_id=compound.id,
                amount=row[amount_col],
                quote_id=q_id,
                supplier="MCule",
                max_lead_time=None,
            )

        self.db.commit()

        return ingredients

    ### SINGLE INSERTION

    def register_compound(
        self,
        *,
        smiles: str,
        bases: list[Compound] | list[int] | None = None,
        tags: None | list = None,
        metadata: None | dict = None,
        return_compound: bool = True,
        commit: bool = True,
        alias: str | None = None,
        return_duplicate: bool = False,
        register_base_if_duplicate: bool = True,
        radical: str = "warning",
        debug: bool = False,
    ) -> Compound:
        """Use a smiles string to add a compound to the database. If it already exists return the compound

        :param smiles: The SMILES string of the compound
        :param bases: A list of :class:`.Compound` objects or IDs that this compound is based on, defaults to ``None``
        :param tags: A list of tags to assign to this compound, defaults to ``None``
        :param metadata: A dictionary of metadata to assign to this compound, defaults to ``None``
        :param return_compound: return the :class:`.Compound` object instead of the integer ID, defaults to ``True``
        :param commit: Commit the changes to the :class:`.Database`, defaults to ``True``
        :param alias: The string alias of this compound, defaults to ``None``
        :param return_duplicate: If ``True`` returns a boolean indicating if this compound previously existed, defaults to ``False``
        :param register_base_if_duplicate: If this compound exists in the :class:`.Database` modify it's ``base`` property, defaults to ``True``
        :param radical: Define the behaviour for dealing with radical atoms in the SMILES. See :class:`.sanitise_smiles`. Defaults to ``'warning'``
        :param debug: Increase verbosity of output, defaults to ``False``
        :returns: The registered/existing :class:`.Compound` object or its ID (depending on ``return_compound``), and optionally a boolean to indicate duplication see ``return_duplicate``
        """

        assert smiles
        assert isinstance(smiles, str), f"Non-string {smiles=}"

        try:
            smiles = sanitise_smiles(
                smiles, sanitisation_failed="error", radical=radical, verbosity=debug
            )
        except SanitisationError as e:
            logger.error(f"Could not sanitise {smiles=}")
            logger.error(str(e))
            return None
        except AssertionError:
            logger.error(f"Could not sanitise {smiles=}")
            return None

        if bases:
            bases = [b.id if isinstance(b, Compound) else b for b in bases]

        inchikey = inchikey_from_smiles(smiles)

        if debug:
            logger.var("inchikey", inchikey)

        compound_id = self.db.insert_compound(
            smiles=smiles,
            inchikey=inchikey,
            tags=tags,
            metadata=metadata,
            warn_duplicate=False,
            commit=False,
            alias=alias,
        )

        duplicate = not bool(compound_id)

        def _return(compound, duplicate, return_compound, return_duplicate):
            if commit:
                self.db.commit()
            if not return_compound and not isinstance(compound, int):
                compound = compound.id
            if return_duplicate:
                return compound, duplicate
            else:
                return compound

        def check_smiles(compound_id, smiles):
            assert compound_id
            db_smiles = self.db.select_where(
                table="compound", query="compound_smiles", key="id", value=compound_id
            )
            (db_smiles,) = db_smiles
            if db_smiles != smiles:
                logger.warning(
                    f"SMILES changed during compound registration: {smiles} --> {db_smiles}"
                )

        def insert_bases(bases, compound_id):
            bases = [b for b in bases if b is not None] or []
            for base in bases:
                self.db.insert_scaffold(
                    base=base,
                    superstructure=compound_id,
                    warn_duplicate=False,
                    commit=False,
                )

        if return_compound or metadata or alias or tags:
            if not compound_id:
                compound = self.compounds[inchikey]

                check_smiles(compound.id, smiles)

            else:
                compound = self.compounds[compound_id]

            if metadata:
                compound.metadata.update(metadata)

            if alias:
                compound.alias = alias

            if tags:
                for tag in tags:
                    compound.tags.add(tag, commit=False)

            if bases and not (not register_base_if_duplicate and duplicate):
                insert_bases(bases, compound.id)

            return _return(compound, duplicate, return_compound, return_duplicate)

        else:
            if not compound_id:
                compound_id = self.db.get_compound_id(inchikey=inchikey)

                check_smiles(compound_id, smiles)

            if bases and not (not register_base_if_duplicate and duplicate):
                insert_bases(bases, compound_id)

            return _return(compound_id, duplicate, return_compound, return_duplicate)

    def register_reaction(
        self,
        *,
        type: str,
        product: Compound | int,
        reactants: list[Compound | int],
        commit: bool = True,
        product_yield: float = 1.0,
        check_chemistry: bool = False,
    ) -> Reaction:
        """Add a :class:`.Reaction` to the :class:`.Database`. If it already exists return the existing one

        :param type: string indicating the type of reaction
        :param product: The :class:`.Compound` object or ID of the product
        :param reactants: A list of :class:`.Compound` objects or IDs of the reactants
        :param commit: Commit the changes to the :class:`.Database`, defaults to ``True``
        :param product_yield: The fraction of product yielded from this reaction ``0 < product_yield <= 1.0``, defaults to ``1.0``
        :param check_chemistry: check the reaction chemistry, defaults to ``True``
        :returns: The registered :class:`.Reaction`
        """

        ### CHECK REACTION VALIDITY

        if check_chemistry:
            from .chem import check_chemistry, InvalidChemistryError

            if not isinstance(product, Compound):
                product = self.db.get_compound(id=product)

            if not isinstance(reactants, CompoundSet):
                reactants = CompoundSet(self.db, reactants)

            valid = check_chemistry(type, reactants, product)

            if not valid:
                raise InvalidChemistryError(f"{type=}, {reactants.ids=}, {product.id=}")

        ### CHECK FOR DUPLICATES

        if isinstance(product, Compound):
            product = product.id

        reactant_ids = set(v.id if isinstance(v, Compound) else v for v in reactants)

        pairs = self.db.execute(
            f'SELECT reactant_reaction, reactant_compound FROM reactant INNER JOIN reaction ON reactant.reactant_reaction = reaction.reaction_id WHERE reaction_type="{type}" AND reaction_product = {product}'
        ).fetchall()

        if pairs:

            reax_dict = {}
            for reaction_id, reactant_id in pairs:
                if reaction_id not in reax_dict:
                    reax_dict[reaction_id] = set()
                reax_dict[reaction_id].add(reactant_id)

            for reaction_id, reactants in reax_dict.items():
                if reactants == reactant_ids:
                    return self.reactions[reaction_id]

        ### INSERT A NEW REACTION

        reaction_id = self.db.insert_reaction(
            type=type, product=product, commit=commit, product_yield=product_yield
        )

        ### INSERT REACTANTS

        for reactant in reactant_ids:
            self.db.insert_reactant(
                compound=reactant, reaction=reaction_id, commit=commit
            )

        return self.reactions[reaction_id]

    def register_target(
        self,
        name: str,
    ) -> Target:
        """
        Register a new protein :class:`` to the Database

        :param param1: this is a first param
        :param param2: this is a second param
        :returns: this is a description of what is returned
        :raises keyError: raises an exception
        """

        target_id = self.db.insert_target(name=name)

        if not target_id:
            target_id = self.db.get_target_id(name=name)

        return self.db.get_target(id=target_id)

    def register_pose(
        self,
        *,
        compound: Compound | int,
        target: str,
        path: str,
        inchikey: str | None = None,
        alias: str | None = None,
        reference: int | None = None,
        tags: None | list = None,
        metadata: None | dict = None,
        inspirations: None | list[int | Pose] = None,
        return_pose: bool = True,
        energy_score: float | None = None,
        distance_score: float | None = None,
        commit: bool = True,
        overwrite_metadata: bool = True,
        warn_duplicate: bool = True,
        check_RMSD: bool = False,
        RMSD_tolerance: float = 1.0,
        split_PDB: bool = False,
        duplicate_alias: str = "modify",
        resolve_path: bool = True,
    ) -> Pose:
        """Add a :class:`.Pose` to the :class:`.Database`. If it already exists return the pose

        :param compound: The :class:`.Compound` object or ID that this :class:`.Pose` is a conformer of
        :param target: The :class:`.Target` name or ID
        :param path: Path to the :class:`.Pose`'s conformer file (.pdb or .mol)
        :param alias: The string alias of this :class:`.Pose`, defaults to ``None``
        :param reference: Reference :class:`.Pose` to use as the protein conformation for all poses, defaults to ``None``
        :param tags: A list of tags to assign to this compound, defaults to ``None``
        :param metadata: A dictionary of metadata to assign to this compound, defaults to ``None``
        :param inspirations: a list of inspiration :class:`.Pose` objects or ID's, defaults to ``None``
        :param energy_score: assign an energy score to this :class:`.Pose`, defaults to ``None``
        :param distance_score: assign a distance score to this :class:`.Pose`, defaults to ``None``
        :param commit: Commit the changes to the :class:`.Database`, defaults to ``True``
        :param overwrite_metadata: If a duplicate is found, overwrite its metadata, defaults to ``True``
        :param warn_duplicate: Warn if a duplicate :class:`.Pose` exists, defaults to ``True``
        :param check_RMSD: Check the RMSD against existing :class:`.Pose`, defaults to ``False``
        :param RMSD_tolerance: Tolerance for ``check_RMSD`` in Angstrom, defaults to ``1.0``
        :param split_PDB: Register a :class:`.Pose` for every ligand residue in the PDB, defaults to ``False``
        :param duplicate_alias: In the case of a duplicate, define the behaviour for the ``alias`` property, defaults to ``'modify'`` which appends ``_copy`` to the alias. Set to ``error`` to raise an Exception.
        :param resolve_path: Resolve to an absoltue path, default = True.
        :returns: The registered/existing :class:`.Pose` object or its ID (depending on ``return_compound``)
        """

        assert duplicate_alias in ["error", "modify"]

        from molparse import parse

        if split_PDB:

            sys = parse(path, verbosity=False, alternative_site_warnings=False)

            lig_residues = []

            for res in sys.ligand_residues:
                lig_residues += res.split_by_site()

            if len(lig_residues) > 1:

                assert not energy_score
                assert not distance_score

                logger.warning(f"Splitting ligands in PDB: {path}")

                results = []
                for i, res in enumerate(lig_residues):
                    file = str(path).replace(".pdb", f"_hippo_{i}.pdb")

                    split_sys = sys.protein_system

                    for atom in res.atoms:
                        split_sys.add_atom(atom)

                    logger.writing(file)
                    split_sys.write(file, verbosity=False)

                    result = self.register_pose(
                        compound=compound,
                        target=target,
                        path=file,
                        inchikey=inchikey,
                        alias=alias,
                        reference=reference,
                        tags=tags,
                        metadata=metadata,
                        inspirations=inspirations,
                        return_pose=return_pose,
                        commit=commit,
                        overwrite_metadata=overwrite_metadata,
                        warn_duplicate=warn_duplicate,
                        check_RMSD=check_RMSD,
                        RMSD_tolerance=RMSD_tolerance,
                        split_PDB=False,
                    )

                    results.append(result)

                return results

        if isinstance(compound, int):
            compound_id = compound
        else:
            compound_id = compound.id

        if check_RMSD:

            # check if the compound has existing poses
            other_pose_ids = self.db.select_id_where(
                table="pose",
                key="compound",
                value=compound_id,
                none="quiet",
                multiple=True,
            )

            if other_pose_ids:
                other_poses = PoseSet(self.db, [i for i, in other_pose_ids])

                from molparse.rdkit import draw_mols, draw_flat
                from rdkit.Chem import MolFromMolFile
                from numpy.linalg import norm
                from numpy import array

                mol = MolFromMolFile(str(path.resolve()))

                c1 = mol.GetConformer()
                atoms1 = [a for a in mol.GetAtoms()]
                symbols1 = [a.GetSymbol() for a in atoms1]
                positions1 = [c1.GetAtomPosition(i) for i, _ in enumerate(atoms1)]

                for pose in other_poses:
                    c2 = pose.mol.GetConformer()
                    atoms2 = [a for a in pose.mol.GetAtoms()]
                    symbols2 = [a.GetSymbol() for a in atoms2]
                    positions2 = [c2.GetAtomPosition(i) for i, _ in enumerate(atoms2)]

                    for s1, p1 in zip(symbols1, positions1):
                        for s2, p2 in zip(symbols2, positions2):
                            if s2 != s1:
                                continue
                            if norm(array(p2 - p1)) <= RMSD_tolerance:
                                # this atom (1) is within tolerance
                                break
                        else:
                            # this atom (1) is outside of tolerance
                            break
                    else:
                        # all atoms within tolerance --> too similar
                        logger.warning(f"Found similar {pose=}")
                        if return_pose:
                            return pose
                        else:
                            return pose.id

        pose_data = dict(
            compound=compound,
            inchikey=inchikey,
            alias=alias,
            target=target,
            path=path,
            tags=tags,
            metadata=metadata,
            reference=reference,
            warn_duplicate=warn_duplicate,
            commit=commit,
            energy_score=energy_score,
            distance_score=distance_score,
        )

        pose_id = self.db.insert_pose(**pose_data, resolve_path=resolve_path)

        # if no pose_id then there must be a duplicate
        if not pose_id:

            # constraint failed
            if isinstance(path, Path):
                path = path.resolve()

            # try getting by path
            result = self.db.select_where(
                table="pose", query="pose_id", key="path", value=str(path), none="quiet"
            )

            # try getting by alias
            if not result:
                result = self.db.select_where(
                    table="pose", query="pose_id", key="alias", value=alias
                )

                if result and duplicate_alias == "error":
                    raise Exception("could not register pose with existing alias")

                elif result and duplicate_alias == "modify":

                    new_alias = alias + "_copy"

                    logger.warning(f"Modifying alias={alias} --> {new_alias}")

                    pose_data["alias"] = new_alias
                    pose_id = self.db.insert_pose(**pose_data)

                else:
                    (pose_id,) = result
            else:
                (pose_id,) = result

            assert pose_id

        if not pose_id:
            logger.var("compound", compound)
            logger.var("inchikey", inchikey)
            logger.var("alias", alias)
            logger.var("target", target)
            logger.var("path", path)
            logger.var("reference", reference)
            logger.var("tags", tags)
            logger.debug(f"{metadata=}")
            logger.debug(f"{inspirations=}")

            raise Exception

        if return_pose or (metadata and not overwrite_metadata):
            pose = self.poses[pose_id]

            if metadata:
                pose.metadata.update(metadata)

        else:
            pose = pose_id

        if overwrite_metadata:
            self.db.insert_metadata(
                table="pose", id=compound_id, payload=metadata, commit=commit
            )

        inspirations = inspirations or []
        for inspiration in inspirations:
            self.db.insert_inspiration(
                original=inspiration,
                derivative=pose,
                warn_duplicate=False,
                commit=commit,
            )

        return pose

    def register_route(
        self,
        *,
        route,
        commit=True,
    ) -> int:
        """
        Insert a :class:`.Route` (single-product :class:`.Recipe`) into the :class:`.Database`.

        :param route: The :class:`.Route` object to be registered
        :param commit: Commit the changes to the :class:`.Database`, defaults to ``True``
        :returns: The :class:`.Route` ID
        """

        # register the route
        route_id = self.db.insert_route(product_id=route.product.id, commit=False)

        assert route_id

        # reactions
        for ref in route.reactions.ids:
            self.db.insert_component(
                component_type=1, ref=ref, route=route_id, commit=False
            )

        # reactants
        for ref in route.reactants.compound_ids:
            self.db.insert_component(
                component_type=2, ref=ref, route=route_id, commit=False
            )

        # intermediates
        for ref in route.intermediates.compound_ids:
            self.db.insert_component(
                component_type=3, ref=ref, route=route_id, commit=False
            )

        if commit:
            self.db.commit()

        return route_id

    ### QUOTING

    def quote_compounds(
        self,
        quoter: Quoter,
        compounds: CompoundSet,
    ) -> None:
        """Get batch quotes using the hippo.Quoter object supplied and add the quotes to the database

        :param quoter: The :class:`.Quoter` object to use
        :param compounds: A :class:`.CompoundSet` containing the compounds to be quoted
        """

        logger.header(
            f"Getting {quoter.supplier} quotes for {len(compounds)} compounds"
        )
        batch_size = quoter.batch_size
        for i in range(0, len(compounds), batch_size):
            logger.debug(f"batch {i%batch_size}")
            batch = compounds[i : i + batch_size]
            quoter.get_batch_quote(batch)

    def quote_reactants(
        self,
        quoter: Quoter,
        *,
        unquoted_only: bool = False,
    ) -> None:
        """Get batch quotes for all reactants in the database

        :param quoter: The :class:`.Quoter` object to use
        :param unquoted_only: Only request quotes for unquoted compouds, defaults to ``False``
        """

        if unquoted_only:
            compounds = self.reactants.get_unquoted(supplier=quoter.supplier)
        else:
            compounds = self.reactants

        self.quote_compounds(quoter=quoter, compounds=compounds)

    def quote_intermediates(
        self,
        quoter: Quoter,
    ) -> None:
        """Get batch quotes for all reactants in the database

        :param quoter: The :class:`.Quoter` object to use
        :param unquoted_only: Only request quotes for unquoted compouds, defaults to ``False``
        """

        self.quote_compounds(quoter=quoter, compounds=self.intermediates)

    ### PLOTTING

    def plot_tag_statistics(self, *args, **kwargs) -> "plotly.graph_objects.Figure":
        """Plot an overview of the number of compounds and poses for each tag, see :func:`hippo.plotting.plot_tag_statistics`"""

        if not self.num_tags:
            logger.error("No tagged compounds or poses")
            return
        from .plotting import plot_tag_statistics

        return plot_tag_statistics(self, *args, **kwargs)

    def plot_compound_property(self, prop, **kwargs) -> "plotly.graph_objects.Figure":
        """Plot an arbitrary compound property across the whole dataset, see :func:`hippo.plotting.plot_compound_property`"""
        from .plotting import plot_compound_property

        return plot_compound_property(self, prop, **kwargs)

    def plot_pose_property(self, prop, **kwargs) -> "plotly.graph_objects.Figure":
        """Plot an arbitrary pose property across the whole dataset, see :func:`hippo.plotting.plot_pose_property`"""
        from .plotting import plot_pose_property

        return plot_pose_property(self, prop, **kwargs)

    def plot_interaction_punchcard(
        self, poses=None, subtitle=None, opacity=1.0, **kwargs
    ) -> "plotly.graph_objects.Figure":
        """Plot an interaction punchcard for a set of poses, see :func:`hippo.plotting.plot_interaction_punchcard`"""
        from .plotting import plot_interaction_punchcard

        return plot_interaction_punchcard(
            self, poses=poses, subtitle=subtitle, opacity=opacity, **kwargs
        )

    def plot_residue_interactions(
        self, poses, residue_number, **kwargs
    ) -> "plotly.graph_objects.Figure":
        """Plot an interaction punchcard for a set of poses, see :func:`hippo.plotting.plot_residue_interactions`"""
        from .plotting import plot_residue_interactions

        return plot_residue_interactions(
            self, poses=poses, residue_number=residue_number, **kwargs
        )

    def plot_compound_availability(
        self, compounds=None, **kwargs
    ) -> "plotly.graph_objects.Figure":
        """Plot a bar chart of compound availability by supplier/catalogue, see :func:`hippo.plotting.plot_compound_availability`"""
        from .plotting import plot_compound_availability

        return plot_compound_availability(self, compounds=compounds, **kwargs)

    def plot_compound_price(
        self,
        min_amount,
        compounds=None,
        plot_lead_time=False,
        style="histogram",
        **kwargs,
    ) -> "plotly.graph_objects.Figure":
        """Plot a bar chart of minimum compound price for a given minimum amount, see :func:`hippo.plotting.plot_compound_price`"""
        from .plotting import plot_compound_price

        return plot_compound_price(
            self, min_amount=min_amount, compounds=compounds, style=style, **kwargs
        )

    def plot_reaction_funnel(self, **kwargs) -> "plotly.graph_objects.Figure":
        """Plot a funnel chart of the reactants, intermediates, and products across the whole dataset, see :func:`hippo.plotting.plot_reaction_funnel`"""
        from .plotting import plot_reaction_funnel

        return plot_reaction_funnel(self, **kwargs)

    def plot_pose_interactions(self, pose: "Pose") -> "plotly.graph_objects.Figure":
        """3d figure showing the interactions between a :class:`.Pose` and the protein. see :func:`hippo.plotting.plot_pose_interactions`"""
        from .plotting import plot_pose_interactions

        return plot_pose_interactions(self, **kwargs)

    ### OTHER

    def summary(self) -> None:
        """Print a text summary of this HIPPO"""
        logger.header(self)
        logger.var("db_path", self.db_path)
        logger.var("#compounds", self.num_compounds)
        logger.var("#poses", self.num_poses)
        logger.var("#reactions", self.num_reactions)
        logger.var("#tags", self.num_tags)
        logger.var("tags", self.tags.unique)
        # logger.var('#products', len(self.products))

    ### DUNDERS

    def __repr__(self) -> str:
        """Returns a command line representation of this HIPPO"""
        return f'HIPPO("{self.name}")'


class InvalidRowError(Exception): ...
