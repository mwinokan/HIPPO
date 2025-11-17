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

from .db import Database
from pathlib import Path

from .tools import inchikey_from_smiles, sanitise_smiles, SanitisationError

import mcol
import mrich
from mrich import print

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
    :param copy_from: optionally initialise this animal by copying the :class:`.Database` at this given path, defaults to None
    :returns: :class:`.HIPPO` object
    """

    def __init__(
        self,
        name: str,
        db_path: str | Path,
        copy_from: str | Path | None = None,
        overwrite_existing: bool = False,
        update_legacy: bool = False,
    ):

        mrich.bold("Creating HIPPO animal")

        self._name = name

        mrich.var("name", name, color="arg")

        if not isinstance(db_path, Path):
            db_path = Path(db_path)

        mrich.var("db_path", db_path, color="file")

        self._db_path = db_path

        if copy_from:
            self._db = Database.copy_from(
                source=copy_from,
                destination=self.db_path,
                animal=self,
                update_legacy=update_legacy,
                overwrite_existing=overwrite_existing,
            )
        else:
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
        self._scaffolds = None
        self._elabs = None

        mrich.success("Initialised animal", f"[var_name]{self}")

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
    def scaffolds(self) -> CompoundSet:
        """Returns compounds that are the basis for one or more elaborations"""
        if (
            self._scaffolds is None
            or self._scaffolds["total_changes"] != self.db.total_changes
        ):
            self._scaffolds = dict(
                set=self.compounds.scaffolds, total_changes=self.db.total_changes
            )
        return self._scaffolds["set"]

    @property
    def num_elabs(self) -> int:
        """Number of compounds that are an elaboration of an existing scaffold"""
        return len(self.elabs)

    @property
    def num_scaffolds(self) -> int:
        """Number of compounds that are the basis for elaborations"""
        return len(self.scaffolds)

    ### BULK INSERTION

    def add_hits(
        self,
        target_name: str,
        metadata_csv: str | Path,
        aligned_directory: str | Path,
        tags: list | None = None,
        skip: list | None = None,
        debug: bool = False,
        load_pose_mols: bool = False,
    ) -> pd.DataFrame:
        """Load in crystallographic hits from a Fragalysis download or XChemAlign alignment.

        For a Fragalysis download `aligned_directory` and `metadata_csv` should point to the `aligned_files` and `metadata.csv` at the root of the extracted download.
        For an XChemAlign dataset the `aligned_directory` should point to the `aligned_files`.

        :param target_name: Name of this protein :class:`.Target`
        :param metadata_csv: Path to the metadata.csv from the Fragalysis download
        :param aligned_directory: Path to the aligned_files directory from the Fragalysis download
        :param skip: optional list of observation names to skip
        :param debug: bool:  (Default value = False)
        :returns: a DataFrame of metadata

        """

        import re
        from enum import Enum
        import molparse as mp
        from rdkit.Chem import PandasTools
        from .tools import remove_other_ligands

        ### Process arguments

        assert aligned_directory, "aligned_directory must be provided"

        skip = skip or []
        tags = tags or ["hits"]

        if not isinstance(aligned_directory, Path):
            aligned_directory = Path(aligned_directory)

        mrich.var("aligned_directory", aligned_directory)

        ### Register Target

        target = self.register_target(name=target_name)

        ### Determine data format

        class DataFormat(Enum):
            Fragalysis_v2 = 1
            XChemAlign_v2 = 2
            XChemAlign_v3 = 3

            def __str__(self):
                return self.name

        subdirs = list(aligned_directory.glob("*[0-9][0-9][0-9][0-9]*"))

        SUBDIR_PATTERN_FRAGALYSIS = re.compile(r"^.*\d{4}[a-z]$")
        SUBDIR_PATTERN_XCA = re.compile(r"^.*-.\d{4}$")

        fragalysis_subdirs_present = any(
            SUBDIR_PATTERN_FRAGALYSIS.match(subdir.name) for subdir in subdirs
        )
        xca_subdirs_present = any(
            SUBDIR_PATTERN_XCA.match(subdir.name) for subdir in subdirs
        )
        assert (
            fragalysis_subdirs_present ^ xca_subdirs_present
        ), "Unexpected mixed data format"

        if fragalysis_subdirs_present:
            data_format = DataFormat.Fragalysis_v2
        else:

            if any(list(subdir.glob("*_artefacts.pdb")) for subdir in subdirs):
                data_format = DataFormat.XChemAlign_v3
            else:
                data_format = DataFormat.XChemAlign_v2

        mrich.var("data_format", data_format)

        ### Counters

        count_directories_tried = 0
        count_compound_registered = 0
        count_poses_registered = 0

        ### Read metadata

        if data_format is DataFormat.Fragalysis_v2:

            assert metadata_csv, "metadata.csv required"

            meta_df = pd.read_csv(metadata_csv)
            curated_tag_cols = [
                c
                for c in meta_df.columns
                if c
                not in [
                    "Code",
                    "Long code",
                    "Compound code",
                    "Smiles",
                    "Downloaded",
                    "Main status",
                    "GOOD count",
                    "MEDIOCRE count",
                    "BAD count",
                    "RefinementResolution",
                ]
                + GENERATED_TAG_COLS
            ]

            mrich.var("curated_tag_cols", curated_tag_cols)

        ### Parse subdirectories

        match data_format:
            case DataFormat.Fragalysis_v2:

                from .fragalysis import parse_observation_longcode

                sdf_pattern = re.compile(r"^.*\d{4}[a-z].sdf$")

                observations = {}

                for path in list(
                    sorted(aligned_directory.glob(f"*[0-9][0-9][0-9][0-9][a-z]"))
                ):

                    name = path.name

                    if name in skip:
                        continue

                    d = dict(
                        name=name,
                        path=path,
                    )

                    ### SDFs

                    sdfs = []

                    for sdf_path in path.glob("*.sdf"):

                        sdf_name = sdf_path.name

                        # fragalysis SDF
                        if sdf_pattern.match(sdf_name):
                            sdfs.append(sdf_path)

                    if not sdfs:
                        mrich.error(name, "has no compatible SDFs", path)
                        continue

                    elif len(sdfs) > 1:
                        mrich.warning(name, "has multiple compatible SDFs", sdfs)

                    d["sdf"] = sdfs[0]

                    ### PDBs

                    pdbs = [
                        p
                        for p in path.glob("*.pdb")
                        if "_ligand" not in p.name
                        and "_apo" not in p.name
                        and "_hippo" not in p.name
                    ]

                    if not len(pdbs) == 1:
                        mrich.error(name, "has invalid PDBs", pdbs)
                        continue

                    d["pdb"] = pdbs[0]

                    observations[name] = d

                    if debug:
                        print(d)

            case _:

                from .xca import parse_observation_longcode

                observations = {}

                match data_format:
                    case DataFormat.XChemAlign_v2:
                        sdf_pattern = re.compile(
                            r"^.*-.\d{4}_._\d*_\d_.*-.\d{4}\+.\+\d*\+\d_ligand\.sdf$"
                        )
                    case DataFormat.XChemAlign_v3:
                        sdf_pattern = re.compile(
                            r"^.*-.\d{4}_._\d*_._\d_.*-.\d{4}\+.\+\d*\+.\+\d_ligand\.sdf$"
                        )

                for path in list(
                    sorted(aligned_directory.glob(f"*[0-9][0-9][0-9][0-9]"))
                ):

                    name = path.name

                    if name in skip:
                        continue

                    ### Group by SDF

                    sdfs = []

                    for sdf_path in sorted(path.glob("*.sdf")):

                        sdf_name = sdf_path.name

                        if sdf_pattern.match(sdf_name):
                            sdfs.append(sdf_path)

                    if not sdfs:
                        mrich.error(name, "has no compatible SDFs", path)
                        continue

                    for i, sdf in enumerate(sdfs):

                        subname = name + chr(ord("a") + i)

                        d = dict(
                            name=subname,
                            path=path,
                            sdf=sdf,
                        )

                        pdb = path / sdf.name.replace("_ligand.sdf", ".pdb")

                        if not pdb.exists():
                            mrich.error(name, "is missing PDB", pdb)
                            continue

                        d["pdb"] = pdb

                    observations[name] = d

        mrich.var("#valid observations", len(observations))

        n_poses = self.num_poses

        for observation_dict in mrich.track(
            observations.values(), prefix="Adding hits..."
        ):

            path = observation_dict["path"]
            name = observation_dict["name"]
            sdf = observation_dict["sdf"]
            pdb = observation_dict["pdb"]

            if debug:
                mrich.debug("Processing", path)

            count_directories_tried += 1

            # load the SDF
            df = PandasTools.LoadSDF(
                str(sdf), molColName="ROMol", idName="ID", strictParsing=True
            )

            # extract fields
            longcode = df.ID[0]
            mol = df.ROMol[0]

            match data_format:
                case DataFormat.Fragalysis_v2:
                    obs_dict = parse_observation_longcode(longcode)
                case DataFormat.XChemAlign_v2:
                    obs_dict = parse_observation_longcode(longcode)

            if debug:
                mrich.debug(name, longcode)

            # parse the PDB file
            if debug:
                mrich.reading(pdb)
            sys = mp.parse(pdb, verbosity=0)

            # create the single ligand bound pdb
            lig_residues = sys.residues["LIG"]
            if len(lig_residues) > 1 or any(
                r.contains_alternative_sites for r in lig_residues
            ):
                sys = remove_other_ligands(
                    sys, obs_dict["residue_number"], obs_dict["chain"]
                )
                sys.prune_alternative_sites("A", verbosity=0)
                pose_path = str(pdb.resolve()).replace(".pdb", "_hippo.pdb")
                mp.write(pose_path, sys, shift_name=True, verbosity=debug)
            else:
                pose_path = str(pdb.resolve())

            # smiles
            smiles = mp.rdkit.mol_to_smiles(mol)
            smiles = sanitise_smiles(smiles, verbosity=debug)

            # create the molecule / pose
            compound_id = self.db.insert_compound(
                smiles=smiles,
                tags=tags,
                warn_duplicate=debug,
                commit=False,
            )

            if not compound_id:

                inchikey = inchikey_from_smiles(smiles)
                compound = self.compounds[inchikey]

                if not compound:
                    mrich.error(
                        "Compound exists in database but could not be found by inchikey"
                    )
                    mrich.var("smiles", smiles)
                    mrich.var("inchikey", inchikey)
                    mrich.var("observation_shortname", name)
                    raise Exception

            else:
                count_compound_registered += 1
                compound = self.compounds[compound_id]

            # metadata

            match data_format:
                case DataFormat.Fragalysis_v2:

                    meta_row = meta_df[meta_df["Code"] == name]
                    if not len(meta_row):
                        assert longcode
                        meta_row = meta_df[meta_df["Long code"] == longcode]

                    assert len(meta_row)

                    metadata = {"fragalysis_longcode": meta_row["Long code"].values[0]}

                    for tag in GENERATED_TAG_COLS:
                        if tag in meta_row.columns:
                            metadata[tag] = meta_row[tag].values[0]

                    pose_tags = set(tags)

                    for tag in curated_tag_cols:
                        if meta_row[tag].values[0]:
                            pose_tags.add(tag)

                case DataFormat.XChemAlign_v2:
                    metadata = {"xca_longcode": longcode}
                    pose_tags = set(tags)

            pose = self.register_pose(
                compound=compound,
                alias=name,
                target=target.id,
                path=pose_path,
                tags=pose_tags,
                metadata=metadata,
                duplicate_alias="skip",
            )

            if load_pose_mols:
                try:
                    pose.mol
                except Exception as e:
                    mrich.error("Could not load molecule", pose)
                    mrich.error(e)

        mrich.var("#directories parsed", count_directories_tried)
        mrich.var("#compounds registered", count_compound_registered)
        mrich.var("#poses registered", self.num_poses - n_poses)

    def load_sdf(
        self,
        *,
        target: str,
        path: str | Path,
        reference: int | Pose | None = None,
        inspirations: list[int] | PoseSet | None = None,
        compound_tags: None | list[str] = None,
        pose_tags: None | list[str] = None,
        mol_col: str = "ROMol",
        name_col: str | None = "ID",
        inspiration_col: str | None = "ref_mols",
        reference_col: str = "ref_pdb",
        energy_score_col: str = "energy_score",
        distance_score_col: str = "distance_score",
        inspiration_map: None | dict = None,
        convert_floats: bool = True,
        skip_equal_dict: dict | None = None,
        skip_not_equal_dict: dict | None = None,
    ) -> None:
        """Add posed virtual hits from an SDF into the database.

        :param target: Name of the protein :class:`.Target`
        :param path: Path to the SDF
        :param reference: Optional single reference :class:`.Pose` to use as the protein conformation for all poses, defaults to ``None``
        :param reference_col: Column that contains reference :class:`.Pose` aliases or ID's
        :param compound_tags: List of string Tags to assign to all created compounds, defaults to ``None``
        :param pose_tags: List of string Tags to assign to all created poses, defaults to ``None``
        :param mol_col: Name of the column containing the ``rdkit.ROMol`` ligands, defaults to ``"ROMol"``
        :param name_col: Name of the column containing the ligand name/alias, defaults to ``"ID"``
        :param inspirations: Optional single set of inspirations :class:`.PoseSet` object or list of IDs to assign as inspirations to all inserted poses, defaults to ``None``
        :param inspiration_col: Name of the column containing the list of inspiration :class:`.Pose` names or ID's, defaults to ``"ref_mols"``
        :param inspiration_map: Optional dictionary or callable mapping between inspiration strings found in ``inspiration_col`` and :class:`.Pose` ids
        :param energy_score_col: Name of the column containing the list of energy scores ``"energy_score"``
        :param distance_score_col: Name of the column containing the list of distance scores, defaults to ``"distance_score"``
        :param convert_floats: Try to convert all values to ``float``, defaults to ``True``
        :param skip_equal_dict: Skip rows where ``any(row[key] == value for key, value in skip_equal_dict.items())``, defaults to ``None``
        :param skip_not_equal_dict: Skip rows where ``any(row[key] != value for key, value in skip_not_equal_dict.items())``, defaults to ``None``

        All non-name columns are added to the Pose metadata.
        N.B. separate .mol files are not created. The molecule binary will only be stored in the .sqlite file and fake paths are added to the database.
        """

        if not isinstance(path, Path):
            path = Path(path)

        skip_equal_dict = skip_equal_dict or {}
        skip_not_equal_dict = skip_not_equal_dict or {}

        mrich.debug(f"{path=}")

        compound_tags = compound_tags or []
        pose_tags = pose_tags or []

        from rdkit.Chem import PandasTools, MolToMolFile, MolFromMolFile
        from molparse.rdkit import mol_to_smiles, mol_to_pdb_block
        from numpy import isnan
        from pandas import read_pickle
        from tempfile import NamedTemporaryFile

        if path.name.endswith(".sdf"):
            df = PandasTools.LoadSDF(str(path.resolve()))
        else:
            df = read_pickle(path)

        df_columns = list(df.columns)

        target = self.register_target(target)

        assert mol_col in df_columns, f"{mol_col=} not in {df_columns}"

        if name_col:
            assert name_col in df_columns, f"{name_col=} not in {df_columns}"

        if inspiration_col and not inspirations:
            assert (
                inspiration_col in df_columns
            ), f"{inspiration_col=} not in {df_columns}"

        if not reference and reference_col:
            assert reference_col in df_columns, f"{reference_col=} not in {df_columns}"

        output_directory = str(path.name).removesuffix(".sdf")
        output_directory = Path(output_directory)
        if not output_directory.exists:
            mrich.writing(f"Creating output directory {output_directory}")
            os.system(f"mkdir -p {output_directory}")

        n_poses = self.num_poses
        n_comps = self.num_compounds

        ### FILTER DATAFRAME

        mrich.var("SDF entries (pre-filter)", len(df))

        df = df[df["ID"] != "ver_1.2"]

        for k, v in skip_equal_dict.items():
            df = df[df[k] == v]

        for k, v in skip_not_equal_dict.items():
            df = df[df[k] != v]

        mrich.var("SDF entries (post-filter)", len(df))

        ### COMPOUND REGISTRATION

        if "smiles" not in df.columns:
            df["smiles"] = df[mol_col].apply(mol_to_smiles)
        smiles = list(set(df["smiles"].values))
        mrich.debug("#smiles", len(smiles))
        mrich.debug("Registering compounds...")
        pairs = self.register_compounds(smiles=smiles, sanitisation_verbosity=False)
        smiles_lookup = {s1: i for s1, (i, s2) in zip(smiles, pairs)}
        inchi_lookup = self.db.get_compound_inchikey_id_dict(
            inchikeys=smiles_lookup.values()
        )

        df["inchikey"] = df["smiles"].apply(lambda x: smiles_lookup.get(x))
        df["compound_id"] = df["inchikey"].apply(lambda x: inchi_lookup.get(x))
        df["compound_id"] = df["compound_id"].fillna(0).astype(int)

        if n := len(df[df["compound_id"].isna()]):
            mrich.error(n, "invalid compound rows")

        cset = self.compounds[set(i for i in df["compound_id"].values if i)]
        for tag in compound_tags:
            cset.add_tag(tag)

        ### POSE REGISTRATION

        if not inspiration_map:
            inspiration_map = self.db.get_pose_alias_id_dict()

        # dicts: (alias, compound, target, path, metadata, inspirations, tags, reference,)
        data = []

        for i, row in mrich.track(df.iterrows(), prefix="Reading SDF rows..."):

            if name_col:
                name = row[name_col].strip() or f"pose_{i}"
                alias = name
            else:
                name = f"pose_{i}"
                alias = None

            mol = row[mol_col]
            inchikey = row["inchikey"]
            smiles = row["smiles"]
            compound_id = row["compound_id"]
            if not compound_id:
                mrich.error("Skipping invalid compound", i)
                continue
            pose_path = (output_directory / f"{name}.fake.mol").resolve()
            energy_score = float(row[energy_score_col])
            distance_score = float(row[distance_score_col])

            # inspirations

            inspiration_list = []

            if isinstance(inspirations, PoseSet):
                inspiration_list = list(inspirations.ids)

            elif inspirations or inspiration_col:

                if inspirations:
                    insp_str = inspirations
                else:
                    insp_str = row[inspiration_col]

                if isinstance(insp_str, str):
                    insp_str = insp_str.removeprefix("[")
                    insp_str = insp_str.removesuffix("]")
                    insp_str = insp_str.replace("'", "")
                    generator = insp_str.split(",")

                elif isinstance(insp_str, float):
                    generator = []

                else:
                    generator = insp_str

                for insp in generator:
                    insp = insp.strip()

                    try:
                        pose_id = int(insp)
                        inspirations.append(pose_id)

                    except ValueError:
                        if (
                            isinstance(inspiration_map, dict)
                            and insp in inspiration_map
                        ):
                            pose_id = inspiration_map[insp]
                            if pose_id:
                                inspiration_list.append(pose_id)
                        elif hasattr(inspiration_map, "__call__"):
                            pose_id = inspiration_map(insp)
                            if pose_id:
                                inspiration_list.append(pose_id)
                        else:
                            mrich.error(
                                f"Could not find inspiration pose with alias={insp}"
                            )
                            continue

            if not reference:
                ref_str = row.get(reference_col)
                if ref_str:
                    try:
                        reference = int(ref_str)
                    except ValueError:
                        reference = inspiration_map[ref_str]
                else:
                    reference = None

            elif isinstance(reference, Pose):
                reference = reference.id

            # metadata
            metadata = {}
            skip = {
                "smiles",
                "inchikey",
                "compound_id",
                inspiration_col,
                name_col,
                mol_col,
                energy_score_col,
                distance_score_col,
                "target_id",
                "reference_id",
                "path",
                "exports",
            }

            for col in df_columns:
                value = row[col]

                if col in skip:
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
                        mrich.warning(f"Skipping metadata from column={col}.")
                    continue

                metadata[col] = value

            data.append(
                dict(
                    alias=alias,
                    compound_id=compound_id,
                    target_id=target.id,
                    path=pose_path,
                    metadata=metadata,
                    inspiration_ids=inspiration_list,
                    reference_id=reference,
                    mol=mol,
                    inchikey=inchikey,
                    smiles=smiles,
                    energy_score=energy_score,
                    distance_score=distance_score,
                )
            )

        ### ACTUALLY DO THE BULK INSERTION

        mrich.debug("Registering poses...")
        ids = self.db.register_poses(data)
        pset = self.poses[ids]
        mrich.debug("Adding tags...")
        for tag in pose_tags:
            pset.add_tag(tag)

        if n := self.num_compounds - n_comps:
            f = mrich.success
        else:
            f = mrich.warning

        f(f"{n} new compounds from {path}")

        if n := self.num_poses - n_poses:
            f = mrich.success
        else:
            f = mrich.warning

        f(f"{n} new poses from {path}")

    def add_syndirella_scaffolds(
        self,
        output_directory: str | Path,
        *,
        pattern: str = "*-*-?-scaffold-check/scaffold-*",
        tags: None | list[str] = None,
        target: int | str = 1,
        debug: bool = False,
    ) -> None:
        """
        Load Poses from Syndirella "scaffold-check" outputs

        :param df_path: Path to the pickled DataFrame or SDF.
        :param tags: list of tags to assign to compounds and poses, defaults to ``None``
        :param target: :class:`.Target` ID or name
        :param pattern: UNIX pattern by which to search for subdirectories
        :param debug: Increase verbosity of output, defaults to ``False``
        :returns: None
        """

        import json

        output_directory = Path(output_directory)

        n_poses = self.num_poses

        mrich.warning("Not setting inspirations and references")

        for subdir in mrich.track(
            list(output_directory.glob(pattern)), prefix="Loading scaffolds..."
        ):

            inchikey = subdir.parent.name.replace("-scaffold-check", "")

            compound = self.compounds[inchikey]

            if debug:
                mrich.var("subdir", subdir)
                mrich.var("inchikey", inchikey)
                mrich.var("compound", compound)

            name = subdir.name

            mol_file = subdir / f"{name}.minimised.mol"
            if not mol_file.exists():
                continue

            json_file = subdir / f"{name}.minimised.json"
            if not json_file.exists():
                continue

            metadata = json.load(open(json_file, "rt"))

            if debug:
                mrich.print(metadata)

            energy_score = (
                metadata["Energy"]["bound"]["total_score"]
                - metadata["Energy"]["unbound"]["total_score"]
            )
            distance_score = metadata["mRMSD"]

            tags = tags or ["Syndirella scaffold"]

            self.register_pose(
                path=mol_file,
                compound=compound,
                target=target,
                tags=tags,
                return_pose=False,
            )

        n_poses = self.num_poses - n_poses

        if n_poses:
            mrich.success(f"Added {n_poses} scaffold Poses")
        else:
            mrich.warning(f"Added {n_poses} scaffold Poses")

    def add_syndirella_elabs(
        self,
        df_path: str | Path,
        max_energy_score: float | None = 0.0,
        max_distance_score: float | None = 2.0,
        require_intra_geometry_pass: bool = True,
        reject_flags: list[str] | None = None,
        register_reactions: bool = True,
        dry_run: bool = False,
        scaffold_route: "Route | None" = None,
        scaffold_compound: "Compound | None" = None,
        pose_tags: list[str] | None = None,
        product_tags: list[str] | None = None,
    ) -> "pd.DataFrame":
        """
        Load Syndirella elaboration compounds and poses from a pickled DataFrame

        :param df_path: Path to the pickled DataFrame
        :param max_energy_score: Filter out poses with `∆∆G` above this value
        :param max_distance_score: Filter out poses with `comRMSD` above this value
        :param require_intra_geometry_pass: Filter out poses with falsy `intra_geometry_pass` values
        :param reject_flags: Filter out rows flagged with strings from this list (default = ["one_of_multiple_products", "selectivity_issue_contains_reaction_atoms_of_both_reactants"])
        :param scaffold_route: Supply a known single-step route to the scaffold product to use if scaffold placements are missing
        :param scaffold_compound: Supply a :class:`.Compound` for the scaffold product to use if scaffold placements are missing
        :param dry_run: Don't insert new records into the database (for debugging/testing)
        :param pose_tags: Add these tags to all inserted poses, defaults to ["syndirella_product", "syndirella_placed"]
        :param product_tags: Add these tags to all inserted product compounds, defaults to ["syndirella_product"]
        :returns: annotated DataFrame
        """

        reject_flags = reject_flags or [
            "one_of_multiple_products",
            "selectivity_issue_contains_reaction_atoms_of_both_reactants",
        ]

        pose_tags = pose_tags or ["syndirella_product", "syndirella_placed"]
        product_tags = product_tags or ["syndirella_product"]

        from .syndirella import reactions_from_row

        df_path = Path(df_path)
        mrich.h3(df_path.name)

        mrich.reading(df_path)
        df = pd.read_pickle(df_path)

        # work out number of reaction steps
        num_steps = max(
            [int(s.split("_")[0]) for s in df.columns if "_product_smiles" in s]
        )
        mrich.var("num_steps", num_steps)

        # add is_scaffold row
        df["is_scaffold"] = df[f"{num_steps}_product_name"].str.contains("scaffold")

        ###### PREP ######

        # flags

        present_flags = set()
        for step in range(num_steps):
            step += 1

            for flags in set(df[df[f"{step}_flag"].notna()][f"{step}_flag"].to_list()):
                for flag in flags:
                    present_flags.add(flag)

        if present_flags:
            mrich.warning("Flags in DataFrame:", present_flags)

        for flag in reject_flags:
            if flag in present_flags:
                for step in range(num_steps):
                    step += 1
                    matches = df[f"{step}_flag"].apply(
                        lambda x: flag in x if x is not None else False
                    )
                    mrich.print(
                        "Filtering out",
                        len(df[matches]),
                        "rows from step",
                        step,
                        "due to",
                        flag,
                    )
                    df = df[~matches]

        # poses

        n_null_mol = len(df[df["path_to_mol"].isna()])
        if n_null_mol:
            df = df[df["path_to_mol"].notna()]
            mrich.var("#rows skipped due to null path_to_mol", n_null_mol)

        if not len(df):
            mrich.warning("No valid rows")
            return None

        # inspirations
        inspiration_sets = set(tuple(sorted(i)) for i in df["regarded"])
        if len(inspiration_sets) != 1:
            mrich.error("Varying inspirations not supported")
            return df

        (inspiration_set,) = inspiration_sets
        inspirations = self.poses[inspiration_set]
        assert len(inspirations) == len(inspiration_set)

        # reference
        template_paths = set(df["template"].to_list())
        assert len(template_paths) == 1, "Multiple references not supported"
        (template_path,) = template_paths
        template_path = Path(template_path)
        mrich.var("template_path", template_path)
        base_name = template_path.name.removesuffix(".pdb").removesuffix("_apo-desolv")
        reference = self.poses[base_name]
        assert reference, "Could not determine reference structure"
        mrich.var("reference", reference)

        target = reference.target

        # subset of rows
        scaffold_df = df[df["is_scaffold"]]
        elab_df = df[~df["is_scaffold"]]
        mrich.var("#scaffold entries", len(scaffold_df))
        mrich.var("#elab entries", len(elab_df))

        if not len(scaffold_df) and not scaffold_route and not scaffold_compound:
            mrich.error("No valid scaffold rows")
            return None

        elif scaffold_route:

            ### SUPPLEMENT THE SCAFFOLD ROWS FROM KNOWN ROUTE

            assert scaffold_route.num_reactions == 1

            product = scaffold_route.products[0].compound
            reaction = scaffold_route.reactions[0]

            assert len(reaction.reactants) == 2

            scaffold_dict = {
                "scaffold_smiles": product.smiles,
                "1_reaction": reaction.type,
                "1_r1_smiles": reaction.reactants[0].smiles,
                "1_r2_smiles": reaction.reactants[1].smiles,
                "1_product_smiles": product.smiles,
                "1_product_name": "scaffold",
                "1_single_reactant_elab": False,
                "1_num_atom_diff": 0,
                "is_scaffold": True,
            }

            scaffold_df = pd.DataFrame([scaffold_dict])

            df = pd.concat([scaffold_df, df])

            scaffold_df = df[df["is_scaffold"]]
            elab_df = df[~df["is_scaffold"]]

        elif scaffold_compound:

            ### SUPPLEMENT PARTIAL SCAFFOLD ROWS FROM KNOWN PRODUCT

            scaffold_dict = {
                "scaffold_smiles": scaffold_compound.smiles,
                "is_scaffold": True,
            }

            scaffold_df = pd.DataFrame([scaffold_dict])

            df = pd.concat([scaffold_df, df])

            scaffold_df = df[df["is_scaffold"]]
            elab_df = df[~df["is_scaffold"]]

        if dry_run:
            mrich.error("Not registering records (dry_run)")
            return df

        ###### ELABS ######

        # bulk register compounds

        smiles_cols = [
            c for c in df.columns if c.endswith("_smiles") and c != "scaffold_smiles"
        ]

        for smiles_col in smiles_cols:

            inchikey_col = smiles_col.replace("_smiles", "_inchikey")
            compound_id_col = smiles_col.replace("_smiles", "_compound_id")

            unique_smiles = df[smiles_col].dropna().unique()

            mrich.debug(
                f"Registering {len(unique_smiles)} compounds from column: {smiles_col}"
            )

            values = self.register_compounds(
                smiles=unique_smiles,
                radical=False,
                sanitisation_verbosity=False,
            )

            orig_smiles_to_inchikey = {
                orig_smiles: inchikey
                for orig_smiles, (inchikey, new_smiles) in zip(unique_smiles, values)
            }

            df[inchikey_col] = df[smiles_col].apply(
                lambda x: orig_smiles_to_inchikey.get(x)
            )

            # get associated IDs
            compound_inchikey_id_dict = self.db.get_compound_inchikey_id_dict(
                list(orig_smiles_to_inchikey.values())
            )
            df[compound_id_col] = df[inchikey_col].apply(
                lambda x: compound_inchikey_id_dict.get(x)
            )

        # bulk register reactions

        if register_reactions:
            for step in range(num_steps):

                step += 1

                mrich.debug(f"Registering reactions for step {step}")

                reaction_dicts = []

                for reaction_name, r1_id, r2_id, product_id in df[
                    [
                        f"{step}_reaction",
                        f"{step}_r1_compound_id",
                        f"{step}_r2_compound_id",
                        f"{step}_product_compound_id",
                    ]
                ].values:

                    # skip invalid rows
                    if pd.isna(r1_id) or pd.isna(product_id):
                        mrich.warning("Can't insert reactions for missing scaffold")
                        continue

                    # reactant IDs

                    reactant_ids = set()
                    reactant_ids.add(int(r1_id))

                    if not pd.isna(r2_id):
                        reactant_ids.add(int(r2_id))

                    product_id = int(product_id)

                    # registration data

                    reaction_dicts.append(
                        dict(
                            reaction_name=reaction_name,
                            reactant_ids=reactant_ids,
                            product_id=int(product_id),
                        )
                    )

            reaction_ids = self.register_reactions(
                types=[d["reaction_name"] for d in reaction_dicts],
                product_ids=[d["product_id"] for d in reaction_dicts],
                reactant_id_lists=[d["reactant_ids"] for d in reaction_dicts],
            )

        scaffold_df = df[df["is_scaffold"]]
        elab_df = df[~df["is_scaffold"]]

        # tag product compounds:

        product_ids = list(df[f"{num_steps}_product_compound_id"].dropna().unique())
        products = self.compounds[product_ids]
        for tag in product_tags:
            products.add_tag(tag)

        # bulk register scaffold relationships

        for step in range(num_steps):

            step += 1

            for role in ["r1", "r2", "product"]:

                key = f"{step}_{role}_compound_id"

                mrich.debug(f"Registering scaffold relatonships for {key}")

                scaffold_ids = list(scaffold_df[key].dropna().unique())

                if not scaffold_ids:
                    mrich.warning(
                        "Can't insert scaffold relationships due to missing",
                        key,
                        "for all scaffold rows",
                    )
                    continue

                if len(scaffold_ids) > 1:
                    mrich.error("Multiple scaffold row values in", key)
                    return scaffold_df

                scaffold_id = scaffold_ids[0]
                superstructure_ids = [
                    i for i in elab_df[key].unique() if i != scaffold_id
                ]

                sql = """
                INSERT OR IGNORE INTO scaffold(scaffold_base, scaffold_superstructure)
                VALUES(?1, ?2)
                """

                self.db.executemany(sql, [(scaffold_id, i) for i in superstructure_ids])
                self.db.commit()

        # filter poses

        ok = df

        try:
            if require_intra_geometry_pass:
                mrich.var(
                    "#poses !intra_geometry_pass",
                    len(df[df["intra_geometry_pass"] == False]),
                )
                ok = ok[ok["intra_geometry_pass"] == True]

            if max_energy_score is not None:
                mrich.var(
                    f"#poses ∆∆G > {max_energy_score}",
                    len(df[df["∆∆G"] > max_energy_score]),
                )
                ok = ok[ok["∆∆G"] <= max_energy_score]

            if max_distance_score is not None:
                mrich.var(
                    f"#poses comRMSD > {max_distance_score}",
                    len(df[df["comRMSD"] > max_energy_score]),
                )
                ok = ok[ok["comRMSD"] <= max_distance_score]

        except Exception as e:
            mrich.error("Problem filtering dataframe")
            mrich.error(e)
            return df

        mrich.var("#acceptable poses", len(ok))

        if not len(ok):
            mrich.warning("No valid poses")
            return None

        # bulk register poses

        payload = []

        for i, row in ok.iterrows():

            path = Path(row.path_to_mol).resolve()

            if not path.exists():
                mrich.warning("Skipping pose w/ non-exising file:", path)
                continue

            pose_tuple = (
                int(reference.id),
                str(path),
                int(row[f"{num_steps}_product_compound_id"]),
                int(target.id),
                float(row["∆∆G"]),
                float(row["comRMSD"]),
            )

            payload.append(pose_tuple)

        if not payload:
            mrich.warning("No valid poses")
            return None

        mrich.debug(f"Registering {len(payload)} poses...")

        sql = """
        INSERT OR IGNORE INTO pose(
            pose_reference,
            pose_path,
            pose_compound,
            pose_target,
            pose_energy_score,
            pose_distance_score
        )
        VALUES(?1, ?2, ?3, ?4, ?5, ?6)
        """

        n_before = self.num_poses
        self.db.executemany(sql, payload)
        self.db.commit()
        diff = self.num_poses - n_before

        if diff:
            mrich.success("Registered", diff, "new poses")
        else:
            mrich.warning("Registered", diff, "new poses")

        # query relevant poses (also previously registered)
        paths = [t[1] for t in payload]
        str_ids = str(tuple(paths)).replace(",)", ")")
        records = self.db.select_where(
            table="pose", query="pose_id", key=f"pose_path IN {str_ids}", multiple=True
        )
        pose_ids = [i for i, in records]

        # bulk register inspirations

        payload = set()
        for pose_id in pose_ids:
            for inspiration in inspirations.ids:
                payload.add((inspiration, pose_id))

        sql = """
        INSERT OR IGNORE INTO inspiration(inspiration_original, inspiration_derivative)
        VALUES(?1, ?2)
        """

        self.db.executemany(sql, list(payload))
        self.db.commit()

        # if pose_tags:
        poses = self.poses[pose_ids]
        for tag in pose_tags:
            poses.add_tag(tag)

        return df

    def add_syndirella_routes(
        self,
        pickle_path: str | Path,
        CAR_only: bool = True,
        pick_first: bool = True,
        check_chemistry: bool = True,
        register_routes: bool = True,
    ) -> pd.DataFrame:
        """Add routes found from syndirella --just_retro query"""

        from .recipe import Recipe
        from .cset import IngredientSet
        from .rset import ReactionSet
        from .chem import InvalidChemistryError, UnsupportedChemistryError

        df = pd.read_pickle(pickle_path)

        for i, row in mrich.track(df.iterrows(), total=len(df)):

            mrich.set_progress_field("i", i)
            mrich.set_progress_field("n", len(df))

            d = row.to_dict()

            comp = self.compounds(smiles=d["smiles"])

            n_routes = 0
            for key in d:
                if not key.startswith("route"):
                    continue

                if not key.endswith("_names"):
                    continue

                v = d[key]

                if isinstance(v, float) and pd.isna(v):
                    break

                n_routes += 1

            if not n_routes:
                # mrich.warning(comp, "#routes =", n_routes)
                continue

            routes = []
            for j in range(n_routes):

                route_str = f"route{j}"

                route = d[route_str]

                if CAR_only and not d[route_str + "_CAR"]:
                    continue

                reactions = ReactionSet(self.db)
                reactants = IngredientSet(self.db)
                intermediates = IngredientSet(self.db)
                products = IngredientSet(self.db)

                try:
                    for k, reaction in enumerate(route):

                        reaction_type = reaction["name"]

                        product = self.compounds(smiles=reaction["productSmiles"])

                        mrich.print(i, j, k, reaction_type, product)

                        rs = []
                        for reactant_s in reaction["reactantSmiles"]:
                            reactant = self.register_compound(smiles=reactant_s)
                            rs.append(reactant.id)

                        # register the reaction
                        reaction = self.register_reaction(
                            type=reaction_type,
                            product=product,
                            reactants=rs,
                            check_chemistry=check_chemistry,
                        )

                        for r_id in rs:
                            if r_id in reactants:
                                intermediates.add(compound_id=r_id, amount=1)
                            else:
                                reactants.add(compound_id=r_id, amount=1)

                        reactions.add(reaction)
                except InvalidChemistryError:
                    continue
                except UnsupportedChemistryError:
                    mrich.warning("Skipping unsupported chemistry:", reaction_type)
                    continue
                except Exception as e:
                    mrich.error("Uncaught error with row", i, "route", j, "reaction", k)
                    continue

                products.add(product.as_ingredient(amount=1))

                recipe = Recipe(
                    db=self.db,
                    reactions=reactions,
                    reactants=reactants,
                    intermediates=intermediates,
                    products=products,
                )

                if register_routes:
                    route_id = self.register_route(recipe=recipe)
                    mrich.success("registered route", route_id)

                if pick_first:
                    break

        return df

    def add_enamine_quote(
        self,
        path: str | Path,
        *,
        orig_name_col: str = "Customer Code",
        # orig_name_col: str = 'Diamond ID (Molecule Name)',
        price_col: str | None = None,
        fixed_amount: float | None = None,
        fixed_lead_time: float | None = False,
        fixed_purity: float | None = False,
        entry_col: str = "Catalog ID",
        catalogue_col: str = "Collection",
        smiles_col: str = "SMILES",
        amount_col: str = "Amount, mg",
        purity_col: str = "Purity, %",
        lead_time_col: str | None = "Lead time",
        stop_after: None | int = None,
        orig_name_is_hippo_id: bool = False,
        allow_no_catalogue_col: bool = False,
        delete_unavailable: bool = True,
        overwrite_existing_quotes: bool = False,
        supplier_name: str = "Enamine",
    ):
        """
        Load an Enamine quote provided as an excel file

        :param path: Path to the excel file
        :param orig_name_col: Column name of the original alias, defaults to 'Customer Code'
        :param entry_col: Column name of the catalogue ID/entry, defaults to 'Catalog ID'
        :param price_col: Column name of the price, defaults to 'Price, EUR' or 'Price, USD' if present
        :param catalogue_col: Column name of the price, defaults to 'Price, EUR' or 'Price, USD' if present
        :param fixed_amount: Optionally use a fixed amount for all quotes (in mg)
        :param fixed_lead_time: Optionally use a fixed lead time for all quotes (in days)
        :param stop_after: Stop after given number of rows, defaults to ``None``
        :param orig_name_is_hippo_id: Set to ``True`` if ``orig_name_col`` is the original HIPPO :class:``hippo.compound.Compound`` ID, defaults to ``False``
        :param delete_unavailable: Delete existing Enamine database quotes for compounds that are unavailable in the quote being loaded
        :param overwrite_existing_quotes: Delete existing Enamine database quotes for compounds that are available in the quote being loaded
        :returns: An :class:`.IngredientSet` of the quoted molecules
        """

        df = pd.read_excel(path)

        def unexpected_column(key, value) -> str:
            return (
                f"Unexpected Excel format ({key}='{value}') \n\nfirst row:\n{df.loc[0]}"
            )

        if smiles_col not in df.columns:
            smiles_col = smiles_col.lower()

        assert smiles_col in df.columns, unexpected_column("smiles_col", smiles_col)

        if orig_name_col is not None:
            assert orig_name_col in df.columns, unexpected_column(
                "orig_name_col", orig_name_col
            )
        else:
            orig_name_is_hippo_id = False

        assert entry_col in df.columns, unexpected_column("entry_col", entry_col)

        if fixed_purity is False:
            assert purity_col in df.columns, unexpected_column("purity_col", purity_col)

        if fixed_amount is None:
            assert amount_col in df.columns, unexpected_column("amount_col", amount_col)

        if fixed_lead_time is False and lead_time_col is not None:
            assert lead_time_col in df.columns, unexpected_column(
                "lead_time_col", lead_time_col
            )

        if not allow_no_catalogue_col:
            assert catalogue_col in df.columns, unexpected_column(
                "catalogue_col", catalogue_col
            )

        assert (
            "Price, EUR" in df.columns
            or "Price, USD" in df.columns
            or price_col in df.columns
        ), unexpected_column("Price", "")

        if price_col is None:
            price_cols = [c for c in df.columns if c.startswith("Price")]
            assert len(price_cols) == 1
            price_col = price_cols[0]
        currency = price_col.split(", ")[-1]

        ingredients = IngredientSet(self.db)

        if len(df) > 100:
            generator = mrich.track(
                df.iterrows(), prefix="Loading quotes...", total=len(df)
            )
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
                        mrich.error("Compound registration mismatch:")
                        mrich.var("expected_id", expected_id)
                        mrich.var("new_id", compound.id)
                        mrich.var("original_smiles", self.compounds[expected_id].smiles)
                        mrich.var("new_smiles", smiles)

                except ValueError:
                    pass

            if (catalogue := row[catalogue_col]) in [
                "No starting material",
                "Out of stock",
                "Unavailable",
            ]:

                if delete_unavailable:

                    mrich.warning(f"Deleting '{supplier_name}' quotes for", compound)

                    self.db.delete_where(
                        table="quote",
                        key=f"quote_supplier = '{supplier_name}' AND quote_compound = {compound.id}",
                    )

                continue

            if (price := row[price_col]) == 0.0:
                continue

            if fixed_amount is None:
                amount = row[amount_col]
            else:
                amount = fixed_amount

            if fixed_purity is False:
                purity = row[purity_col] / 100
            else:
                purity = fixed_purity

            if fixed_lead_time is False:
                if not isinstance(row[lead_time_col], str):
                    continue
                if "week" in row[lead_time_col]:
                    lead_time = int(row[lead_time_col].split()[0].split("-")[-1]) * 5
                else:
                    raise NotImplementedError
            else:
                lead_time = fixed_lead_time

            if overwrite_existing_quotes:
                self.db.delete_where(
                    table="quote",
                    key=f"quote_supplier = '{supplier_name}' AND quote_compound = {compound.id}",
                )

            quote_data = dict(
                compound=compound,
                supplier=supplier_name,
                catalogue=catalogue,
                entry=row[entry_col],
                amount=amount,
                purity=purity,
                lead_time=lead_time,
                price=price,
                currency=currency,
                smiles=smiles,
            )

            q_id = self.db.insert_quote(**quote_data)

            ingredients.add(
                compound_id=compound.id,
                amount=amount,
                quoted_amount=amount,
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

        for i, row in mrich.track(df.iterrows(), prefix="Loading quotes..."):
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

    def add_soakdb_compounds(
        self,
        path: "str | Path",
        smiles_col: str = "CompoundSMILES",
        alias_col: str = "CompoundCode",
        update_aliases: bool = True,
        soak_count_to_metadata: bool = True,
        sanitisation_verbosity: bool = False,
        stop_after: int | None = None,
    ) -> "CompoundSet":
        """Registers compounds with aliases and metadata from a SoakDB file

        :param path: Path to SoakDB CSV or SQLite file
        :returns: :class:`.CompoundSet` of registered/matched compounds
        """

        from json import dumps

        path = Path(path)

        match ext := path.name.split(".")[-1]:
            case "csv":
                df = pd.read_csv(path)
            case "sqlite":
                raise NotImplementedError
            case _:
                print(ext)
                raise ValueError(
                    f"Could not determine file type from extension, use '.csv' or '.sqlite' {path}"
                )

        smiles_alias_tuples = []

        for i, row in df.iterrows():

            smiles = row[smiles_col]
            alias = row[alias_col]

            if pd.isna(smiles):
                continue

            if pd.isna(alias):
                continue

            smiles_alias_tuples.append((smiles, alias))

            if stop_after and i > stop_after:
                break

        smiles_alias_tuples = set(smiles_alias_tuples)

        old_smiles = [s for s, a in smiles_alias_tuples]

        mrich.debug("Registering compounds...")
        inchikey_new_smiles_tuples = self.register_compounds(
            smiles=old_smiles, sanitisation_verbosity=sanitisation_verbosity
        )

        inchikey_old_smiles_lookup = {
            inchikey: old_s
            for old_s, (inchikey, new_s) in zip(old_smiles, inchikey_new_smiles_tuples)
        }
        alias_lookup = {s: a for s, a in smiles_alias_tuples}
        alias_dicts = [
            dict(compound_inchikey=inchikey, compound_alias=alias_lookup[old_s])
            for old_s, (inchikey, new_s) in zip(old_smiles, inchikey_new_smiles_tuples)
        ]

        if update_aliases:

            sql = """
            UPDATE OR IGNORE compound
            SET compound_alias = :compound_alias
            WHERE compound_inchikey = :compound_inchikey;
            """

            mrich.debug("Updating aliases...")
            self.db.executemany(sql, alias_dicts)

        inchikeys = [d["compound_inchikey"] for d in alias_dicts]

        inchikey_id_lookup = self.db.get_compound_inchikey_id_dict(inchikeys)

        cset = self.compounds[
            [inchikey_id_lookup[d["compound_inchikey"]] for d in alias_dicts]
        ]

        cset.add_tag("soaks")

        metadata_lookup = self.db.get_id_metadata_dict(table="compound", ids=cset.ids)

        if soak_count_to_metadata:

            mrich.debug("Getting soak counts...")
            for inchikey in inchikeys:
                old_s = inchikey_old_smiles_lookup[inchikey]
                c_id = inchikey_id_lookup[inchikey]
                metadata_lookup[c_id]["SoakDB count"] = len(df[df[smiles_col] == old_s])

            sql = """
            UPDATE compound
            SET compound_metadata = ?2
            WHERE compound_id = ?1;
            """

            mrich.debug("Updating metadata...")
            self.db.executemany(
                sql, [(i, dumps(m)) for i, m in metadata_lookup.items()]
            )

        self.db.commit()

        return cset

    ### REGISTRATION

    def register_compound(
        self,
        *,
        smiles: str,
        scaffolds: list[Compound] | list[int] | None = None,
        tags: None | list = None,
        metadata: None | dict = None,
        return_compound: bool = True,
        commit: bool = True,
        alias: str | None = None,
        return_duplicate: bool = False,
        register_scaffold_if_duplicate: bool = True,
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
        :param register_scaffold_if_duplicate: If this compound exists in the :class:`.Database` modify it's ``base`` property, defaults to ``True``
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
            mrich.error(f"Could not sanitise {smiles=}")
            mrich.error(str(e))
            return None
        except AssertionError:
            mrich.error(f"Could not sanitise {smiles=}")
            return None

        if scaffolds:
            scaffolds = [b.id if isinstance(b, Compound) else b for b in scaffolds]

        inchikey = inchikey_from_smiles(smiles)

        if debug:
            mrich.var("inchikey", inchikey)

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
                mrich.warning(
                    f"SMILES changed during compound registration: {smiles} --> {db_smiles}"
                )

        def insert_scaffolds(scaffolds, compound_id):
            scaffolds = [b for b in scaffolds if b is not None] or []
            for scaffold in scaffolds:
                self.db.insert_scaffold(
                    scaffold=scaffold,
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

            if scaffolds and not (not register_scaffold_if_duplicate and duplicate):
                insert_scaffolds(scaffolds, compound.id)

            return _return(compound, duplicate, return_compound, return_duplicate)

        else:
            if not compound_id:

                assert inchikey

                compound_id = self.db.get_compound_id(inchikey=inchikey)

                check_smiles(compound_id, smiles)

            if scaffolds and not (not register_scaffold_if_duplicate and duplicate):
                insert_scaffolds(scaffolds, compound_id)

            return _return(compound_id, duplicate, return_compound, return_duplicate)

    def register_compounds(
        self,
        *,
        smiles: list[str],
        radical: str = "warning",
        sanitisation_verbosity: bool = True,
        debug: bool = False,
    ) -> list[tuple[str, str]]:
        """Insert many compounds at once

        :param smiles: list of smiles strings
        :returns: list of sanitised inchikey and smiles string pairs

        """

        if debug:
            mrich.var("#smiles", len(smiles))

        n_before = self.num_compounds

        values = self.db.register_compounds(
            smiles=smiles,
            radical=radical,
            sanitisation_verbosity=sanitisation_verbosity,
            debug=debug,
        )

        diff = self.num_compounds - n_before

        if diff:
            mrich.success(f"Inserted {diff} new compounds")
        else:
            mrich.warning(f"Inserted {diff} new compounds")

        return values

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
            from .chem import (
                check_chemistry,
                InvalidChemistryError,
                UnsupportedChemistryError,
            )

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
            f"""SELECT reactant_reaction, reactant_compound 
                FROM reactant INNER JOIN reaction 
                ON reactant.reactant_reaction = reaction.reaction_id 
                WHERE reaction_type="{type}" 
                AND reaction_product = {product}"""
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

        assert (
            product_yield > 0 and product_yield <= 1.0
        ), f"{product_yield=} out of range (0,1)"

        reaction_id = self.db.insert_reaction(
            type=type, product=product, commit=commit, product_yield=product_yield
        )

        ### INSERT REACTANTS

        for reactant in reactant_ids:
            self.db.insert_reactant(
                compound=reactant, reaction=reaction_id, commit=commit
            )

        return self.reactions[reaction_id]

    def register_reactions(
        self,
        *,
        types: list[str],
        product_ids: list[list[int]],
        reactant_id_lists: list[list[int]],
    ):
        """Insert many reactions at once

        :param types: list of reaction type strings
        :param reactant_id_lists: list of reactant compound id lists
        :param product_ids: list of product compound ids
        :returns: list of reaction ids
        """

        assert len(types) == len(reactant_id_lists) == len(product_ids)

        # assert not any(not isinstance(t, str) for t in types)
        # assert not any(not isinstance(t, int) for t in product_ids)
        # assert not any(not any(not isinstance(i, int) for i in t) for t in reactant_id_lists)

        types = [str(t) for t in types]
        product_ids = [int(i) for i in product_ids]
        reactant_id_lists = [set(int(i) for i in r) for r in reactant_id_lists]

        n_before = self.num_reactions

        # get possible duplicates
        existing = self.db.get_reaction_map_from_products(product_ids)

        non_duplicates = {}

        existing_count = 0

        for reaction_type, product_id, reactant_ids in zip(
            types, product_ids, reactant_id_lists
        ):

            key = (reaction_type, product_id)

            reactant_ids = set(reactant_ids)

            possible_matches = {k: v for k, v in existing.items() if k == key}

            assert len(possible_matches) < 2

            if possible_matches:
                possible_matches = list(possible_matches.values())[0]

            if any(reactant_ids == v for v in possible_matches.values()):
                existing_count += 1
                continue

            non_duplicates[key] = reactant_ids

        if existing_count:
            mrich.warning("Skipped", existing_count, "existing reactions")

        if not non_duplicates:
            mrich.warning("All reactions are duplicates")
            return None

        # insert reaction records
        sql = """
        INSERT INTO reaction(reaction_type, reaction_product, reaction_product_yield)
        VALUES(?1, ?2, 1)
        RETURNING reaction_id
        """

        payload = list(non_duplicates.keys())

        records = self.db.executemany(sql, payload)
        reaction_ids = [r_id for r_id, in records]
        self.db.commit()

        # insert reactant records
        sql = """
        INSERT OR IGNORE INTO reactant(reactant_amount, reactant_reaction, reactant_compound)
        VALUES(1.0, ?1, ?2)
        """

        payload = []
        for reaction_id, ((reaction_type, product_id), reactant_ids) in zip(
            reaction_ids, non_duplicates.items()
        ):
            for reactant_id in reactant_ids:
                payload.append((reaction_id, reactant_id))

        self.db.executemany(sql, payload)
        self.db.commit()

        diff = self.num_reactions - n_before

        # delete orphaned reactions

        sql = """
        SELECT reaction_id FROM reaction
        LEFT JOIN reactant ON reaction_id = reactant_reaction
        WHERE reactant_compound IS NULL
        """

        records = self.db.execute(sql).fetchall()
        orphaned_str_ids = str(tuple(r for r, in records)).replace(",)", ")")
        self.db.execute(f"DELETE FROM reaction WHERE reaction_id IN {orphaned_str_ids}")

        if diff:
            mrich.success(f"Inserted {diff} new reactions")
        else:
            mrich.warning(f"Inserted {diff} new reactions")

        return reaction_ids

    def register_target(
        self,
        name: str,
        warn_duplicate: bool = True,
    ) -> Target:
        """
        Register a new protein :class:`` to the Database

        :param param1: this is a first param
        :param param2: this is a second param
        :returns: this is a description of what is returned
        :raises keyError: raises an exception
        """

        target_id = self.db.insert_target(name=name, warn_duplicate=warn_duplicate)

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
        load_mol: bool = False,
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
        :param load_mol: Parse the input file and load the ligand rdkit.Chem.Mol
        :returns: The registered/existing :class:`.Pose` object or its ID (depending on ``return_pose``)
        """

        assert duplicate_alias in ["error", "modify", "skip"]

        from molparse import parse

        if split_PDB:

            sys = parse(path, verbosity=False, alternative_site_warnings=False)

            lig_residues = []

            for res in sys.ligand_residues:
                lig_residues += res.split_by_site()

            if len(lig_residues) > 1:

                assert not energy_score
                assert not distance_score

                mrich.warning(f"Splitting ligands in PDB: {path}")

                results = []
                for i, res in enumerate(lig_residues):
                    file = str(path).replace(".pdb", f"_hippo_{i}.pdb")

                    split_sys = sys.protein_system

                    for atom in res.atoms:
                        split_sys.add_atom(atom)

                    mrich.writing(file)
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
                        load_mol=load_mol,
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
                        mrich.warning(f"Found similar {pose=}")
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

                    mrich.warning(f"Modifying alias={alias} --> {new_alias}")

                    pose_data["alias"] = new_alias
                    pose_id = self.db.insert_pose(**pose_data)

                elif result and duplicate_alias == "skip":
                    (pose_id,) = result

                else:
                    (pose_id,) = result
            else:
                (pose_id,) = result

            assert pose_id, (result, pose_id)

        if not pose_id:
            mrich.var("compound", compound)
            mrich.var("inchikey", inchikey)
            mrich.var("alias", alias)
            mrich.var("target", target)
            mrich.var("path", path)
            mrich.var("reference", reference)
            mrich.var("tags", tags)
            mrich.debug(f"{metadata=}")
            mrich.debug(f"{inspirations=}")

            raise Exception

        if return_pose or (metadata and not overwrite_metadata) or load_mol:
            pose = self.poses[pose_id]

            if metadata:
                pose.metadata.update(metadata)

            if load_mol:
                pose.mol

        else:
            pose = pose_id

        if overwrite_metadata:
            self.db.insert_metadata(
                table="pose", id=pose_id, payload=metadata, commit=commit
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
        recipe: "Recipe",
        commit: bool = True,
    ) -> int:
        """
        Insert a single-product :class:`.Recipe` into the :class:`.Database`.

        :param recipe: The :class:`.Recipe` object to be registered
        :param commit: Commit the changes to the :class:`.Database`, defaults to ``True``
        :returns: The :class:`.Route` ID
        """

        assert recipe.num_products == 1

        # register the route
        route_id = self.db.insert_route(product_id=recipe.product.id, commit=False)

        assert route_id

        # reactions
        for ref in recipe.reactions.ids:
            self.db.insert_component(
                component_type=1, ref=ref, route=route_id, commit=False
            )

        # reactants
        for ref, amount in recipe.reactants.id_amount_pairs:
            self.db.insert_component(
                component_type=2, ref=ref, route=route_id, amount=amount, commit=False
            )

        # intermediates
        for ref, amount in recipe.intermediates.id_amount_pairs:
            self.db.insert_component(
                component_type=3, ref=ref, route=route_id, amount=amount, commit=False
            )

        if commit:
            self.db.commit()

        return route_id

    ### QUOTING

    def quote_compounds(
        self,
        ref_animal: "HIPPO",
        compounds: CompoundSet | None = None,
        *,
        debug: bool = False,
    ) -> "CompoundSet,CompoundSet":
        """Transfer quotes from another reference :class:`.HIPPO` animal object (e.g. the one from https://github.com/mwinokan/EnamineCatalogs)

        :param ref_animal: The reference :class:`.HIPPO` animal to fetch quotes from
        :param compounds: A :class:`.CompoundSet` containing the compounds to be quoted
        """

        if compounds is not None:
            inchikeys = compounds.inchikeys

        else:
            inchikeys = self.compounds.inchikeys

        sql = f"""
        SELECT quote_id, quote_smiles, quote_amount, quote_supplier, quote_catalogue, quote_entry, quote_lead_time, quote_price, quote_currency, quote_purity, quote_date, quote_compound FROM quote
        INNER JOIN compound ON quote_compound = compound_id
        WHERE compound_inchikey IN {tuple(inchikeys)}
        """

        with mrich.loading("Querying reference database..."):
            records = ref_animal.db.execute(sql).fetchall()

        quoted_compound_ids = set()
        quote_count = self.db.count("quote")

        for record in mrich.track(
            records, total=len(records), prefix="Inserting quotes"
        ):

            (
                quote_id,
                quote_smiles,
                quote_amount,
                quote_supplier,
                quote_catalogue,
                quote_entry,
                quote_lead_time,
                quote_price,
                quote_currency,
                quote_purity,
                quote_date,
                quote_compound,
            ) = record

            try:
                compound = self.compounds(smiles=quote_smiles)
            except Exception as e:
                mrich.error(e)
                continue

            if debug:
                mrich.debug("Inserting quote for", compound)

            try:
                self.db.insert_quote(
                    compound=compound,
                    supplier=quote_supplier,
                    catalogue=quote_catalogue,
                    entry=quote_entry,
                    amount=quote_amount,
                    price=quote_price,
                    currency=quote_currency,
                    purity=quote_purity,
                    lead_time=quote_lead_time,
                    smiles=quote_smiles,
                    date=quote_date,
                    commit=False,
                )
            except Exception as e:
                mrich.error(e)
                continue

            quoted_compound_ids.add(compound.id)

        self.db.commit()

        quoted_compounds = self.compounds[quoted_compound_ids]

        if compounds:
            unquoted_compounds = compounds - quoted_compounds
        else:
            unquoted_compounds = self.compounds[:] - quoted_compounds

        mrich.var("#new quotes", self.db.count("quote") - quote_count)
        mrich.var("#quoted_compounds", len(quoted_compounds))
        mrich.var("#unquoted_compounds", len(unquoted_compounds))

        return quoted_compounds, unquoted_compounds

    def quote_reactants(
        self,
        ref_animal: "HIPPO",
        *,
        unquoted_only: bool = False,
        supplier: str = "any",
        debug: bool = False,
    ) -> None:
        """Get batch quotes for all reactants in the database

        :param ref_animal: The reference :class:`.HIPPO` animal to fetch quotes from (e.g. the one from https://github.com/mwinokan/EnamineCatalogs)
        :param unquoted_only: Only request quotes for unquoted compouds, defaults to ``False``
        """

        if unquoted_only:
            compounds = self.reactants.get_unquoted(supplier=supplier)
        else:
            compounds = self.reactants

        mrich.var("#compounds", len(compounds))

        self.quote_compounds(ref_animal=ref_animal, compounds=compounds, debug=debug)

    def quote_intermediates(
        self,
        ref_animal: "HIPPO",
    ) -> None:
        """Get batch quotes for all reactants in the database

        :param ref_animal: The reference :class:`.HIPPO` animal to fetch quotes from (e.g. the one from https://github.com/mwinokan/EnamineCatalogs)
        :param unquoted_only: Only request quotes for unquoted compouds, defaults to ``False``
        """

        self.quote_compounds(quoter=quoter, compounds=self.intermediates)

    ### PLOTTING

    def plot_tag_statistics(self, *args, **kwargs) -> "plotly.graph_objects.Figure":
        """Plot an overview of the number of compounds and poses for each tag, see :func:`hippo.plotting.plot_tag_statistics`"""

        if not self.num_tags:
            mrich.error("No tagged compounds or poses")
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

    def plot_interaction_punchcard_by_tags(
        self, tags: dict[str, str] | list[str], **kwargs
    ) -> "plotly.graph_objects.Figure":
        """Plot an interaction punchcard for a set of poses associated to given tags, see :func:`hippo.plotting.plot_interaction_punchcard_by_tags`"""
        from .plotting import plot_interaction_punchcard_by_tags

        return plot_interaction_punchcard_by_tags(self, tags=tags, **kwargs)

    def plot_residue_interactions(
        self, residue_number: int, poses: str | None = None, **kwargs
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

    def plot_compound_availability_venn(
        self, compounds, **kwargs
    ) -> "plotly.graph_objects.Figure":
        """Plot a venn diagram of compound availability by supplier/catalogue, see :func:`hippo.plotting.plot_compound_availability`"""
        from .plotting import plot_compound_availability_venn

        return plot_compound_availability_venn(self, compounds=compounds, **kwargs)

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

    def plot_pose_interactions(
        self, pose: "Pose", **kwargs
    ) -> "plotly.graph_objects.Figure":
        """3d figure showing the interactions between a :class:`.Pose` and the protein. see :func:`hippo.plotting.plot_pose_interactions`"""
        from .plotting import plot_pose_interactions

        return plot_pose_interactions(self, pose=pose, **kwargs)

    def get_scaffold_network(
        self,
        compounds: "CompoundSet | None" = None,
        scaffolds: "CompoundSet | None" = None,
        notebook: bool = True,
        depth: int = 5,
        scaffold_tag: str | None = None,
        exclude_tag: str | None = None,
        physics: bool = True,
        arrows: bool = True,
    ) -> "pyvis.network.Network":
        """Use PyVis to display a network of molecules connected by scaffold relationships in the database"""
        from .pyvis import get_scaffold_network

        return get_scaffold_network(
            self,
            compounds=compounds,
            scaffolds=scaffolds,
            notebook=notebook,
            depth=depth,
            scaffold_tag=scaffold_tag,
            exclude_tag=exclude_tag,
            physics=physics,
            arrows=arrows,
        )

    ### COMPOUND DESIGN

    # def fragmenstein_merge(self,
    #     reference: Pose,
    #     hits: PoseSet,
    #     combination_size: int = 2,
    #     timeout=300,
    #     n_cores: int = 1,
    #     scratch_dir: str = "fragmenstein_scratch",
    #     require_outome: str = "acceptable",
    #     return_df: bool = False,
    #     bulkdock_csv: str = "",
    # ) -> PoseSet:

    #     mrich.var("reference", reference)
    #     mrich.var("combination_size", combination_size)
    #     mrich.var("require_outome", require_outome)
    #     mrich.var("n_cores", n_cores)
    #     mrich.var("timeout", timeout)
    #     mrich.var("scratch_dir", scratch_dir)
    #     mrich.var("protein_path", reference.apo_path)
    #     mrich.var("bulkdock_csv", bulkdock_csv)

    #     from .fstein import setup_wictor_laboratory, pure_merge

    #     lab = setup_wictor_laboratory(
    #         scratch_dir=scratch_dir,
    #         protein_path=reference.apo_path,
    #     )

    #     df = pure_merge(
    #         lab,
    #         hits.mols,
    #         n_cores=n_cores,
    #         timeout=timeout,
    #         combination_size=combination_size,
    #     )

    #     if not len(filtered):
    #         mrich.error(f"No merges")
    #         return None

    #     if require_outome:
    #         filtered = df[df["outcome"] == require_outome]

    #         if not len(filtered):
    #             mrich.error(f"No merges with 'outcome' == {require_outome}")
    #             if return_df:
    #                 return None, df
    #             return None

    #         df = filtered

    #     # register the poses
    #     n = len(self.compounds)
    #     compound_ids = []
    #     # inspirations = []
    #     for i,row in df.iterrows():
    #         smiles = row.smiles
    #         compound = self.register_compound(smiles=smiles, tags=["Fragmenstein pure merge"])
    #         compound_ids.append(compound.id)
    #         # hit_mols = row.hit_names
    #         # print(row)
    #         # raise NotImplementedError

    #     compounds = self.compounds[compound_ids]

    #     d = len(self.compounds) - n

    #     if d:
    #         mrich.success(f"Found and registered {d} new merges")
    #     else:
    #         mrich.success(f"Found and registered {d} new merges")

    #     if return_df:
    #         return compounds, df

    #     return compounds

    ### OTHER

    def summary(self) -> None:
        """Print a text summary of this HIPPO"""
        mrich.header(self)
        mrich.var("db_path", self.db_path)
        mrich.var("#compounds", self.num_compounds)
        mrich.var("#poses", self.num_poses)
        mrich.var("#reactions", self.num_reactions)
        mrich.var("#tags", self.num_tags)
        mrich.var("tags", self.tags.unique)
        # mrich.var('#products', len(self.products))

    def get_by_shorthand(self, key) -> "Compound | Pose | Reaction":
        """Get a :class:`.Compound`, :class:`.Pose`, or :class:`.Reaction` by its ID

        :param key: shortname of the object, e.g. C100 for :class:`.Compound` with id=100
        :returns: :class:`.Compound`, :class:`.Pose`, or :class:`.Reaction` object
        """

        assert isinstance(key, str), f"'HIPPO' object has no attribute '{key}'"
        assert len(key) > 1, f"'HIPPO' object has no attribute '{key}'"

        prefix = key[0]
        index = key[1:]

        if prefix not in "CPR":
            raise AttributeError(f"'HIPPO' object has no attribute '{key}'")

        try:
            index = int(index)
        except ValueError:
            mrich.error(f"Cannot convert {index} to integer")
            return None

        match key[0]:
            case "C":
                return self.compounds[index]
            case "P":
                return self.poses[index]
            case "R":
                return self.reactions[index]

        mrich.error(f"Unsupported {prefix=}")
        return None

    ### DUNDERS

    def __str__(self) -> str:
        """Unformatted string representation of this HIPPO"""
        return f'HIPPO("{self.name}")'

    def __repr__(self) -> str:
        """Returns a command line representation of this HIPPO"""
        return f"{mcol.bold}{mcol.underline}{self}{mcol.clear}"

    def __rich__(self) -> str:
        """Representation for mrich"""
        return f"[bold underline]{self}"

    def __getitem__(self, key: str):
        """Get a :class:`.Compound`, :class:`.Pose`, or :class:`.Reaction` by its ID. See :meth:`.HIPPO.get_by_shorthand`"""
        return self.get_by_shorthand(key)

    def __getattr__(self, key: str):
        """Get a :class:`.Compound`, :class:`.Pose`, or :class:`.Reaction` by its ID. See :meth:`.HIPPO.get_by_shorthand`"""
        return self.get_by_shorthand(key)


GENERATED_TAG_COLS = [
    "ConformerSites alias",
    "CanonSites alias",
    "CrystalformSites alias",
    "Quatassemblies alias",
    "Crystalforms alias",
    "ConformerSites upload name",
    "CanonSites upload name",
    "CrystalformSites upload name",
    "Quatassemblies upload name",
    "Crystalforms upload name",
    "ConformerSites short tag",
    "CanonSites short tag",
    "CrystalformSites short tag",
    "Quatassemblies short tag",
    "Crystalforms short tag",
    "Centroid res",
    "Experiment code",
    "Pose",
]


class InvalidRowError(Exception): ...
