import mrich
from .database import Database
from django.core.exceptions import ObjectDoesNotExist
import molparse as mp
import re
from pathlib import Path
import mcol


class HIPPO:
    """
    HIPPO class
    """

    def __init__(self, name: str, db_path: str, **kwargs) -> None:

        self._name = name
        self._db_path = db_path

        # DB must be initialised before importing any models
        self._db = Database(path=db_path, **kwargs)

        from .compound.compound_table import CompoundTable
        from .pose.pose_table import PoseTable
        from .protein.target_table import TargetTable
        from .chemistry.reaction_table import ReactionTable

        self.compounds = CompoundTable()
        self.poses = PoseTable()
        self.targets = TargetTable()
        self.reactions = ReactionTable()

    ### FACTORIES

    ### PROPERTIES

    @property
    def name(self):
        return self._name

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
    def num_compounds(self) -> int:
        """Total number of Compounds in the Database"""
        return len(self.compounds)

    @property
    def num_poses(self) -> int:
        """Total number of poses in the Database"""
        return len(self.poses)

    @property
    def num_targets(self) -> int:
        """Total number of targets in the Database"""
        return len(self.targets)

    ### DATA LOADING

    def add_hits(
        self,
        target_name: str,
        target_access_string: str | None = None,
        download_directory: "str | Path | None" = None,
        metadata_csv: "str | Path | None" = None,
        aligned_dir: "str | Path | None" = None,
        stack: str = "production",
    ):

        import pandas as pd
        from rdkit.Chem import PandasTools
        from .annotation import Tag, TagType
        from .io.fragalysis import POSE_META_TAG_FIELDS

        if target_access_string:

            from .io.fragalysis import download_target
            from django.conf import settings

            # mrich.print(settings.BASE_DIR)

            assert download_directory, "download_directory required"

            download_directory = Path(download_directory)

            if not download_directory.exists():
                download_directory.mkdir()

            target_dir = download_target(
                target_name=target_name,
                target_access_string=target_access_string,
                destination=download_directory,
                stack=stack,
            )

            metadata_csv = target_dir / "metadata.csv"
            aligned_dir = target_dir / "aligned_files"

        else:

            assert metadata_csv
            assert aligned_dir
            metadata_csv = Path(metadata_csv)
            aligned_dir = Path(aligned_dir)

        assert metadata_csv.exists()
        assert aligned_dir.exists()

        meta_df = pd.read_csv(metadata_csv)

        target = self.register_target(name=target_name)

        shortcodes = set(meta_df["Code"].values)

        poses = []

        for subdir in aligned_dir.iterdir():

            if subdir.name not in shortcodes:
                mrich.debug("Skipping", subdir)

            alias = subdir.name

            meta = meta_df[meta_df["Code"] == alias].iloc[0].to_dict()

            # files = list(subdir.iterdir())

            combi_sdf = subdir / f"{alias}.sdf"
            complex_pdb = subdir / f"{alias}.pdb"
            apo_pdb = subdir / f"{alias}_apo-desolv.pdb"

            mol_df = PandasTools.LoadSDF(
                str(combi_sdf),
                molColName="ROMol",
                idName="ID",
                strictParsing=True,
            )

            assert len(mol_df) == 1, "Combi-soaks not yet supported"

            mol = mol_df.ROMol[0]

            pose = self.register_pose(
                smiles=meta["Smiles"],
                mol=mol,
                complex_file=complex_pdb,
                protein_file=apo_pdb,
                mol_file=combi_sdf,
                alias=alias,
                origin="EXPERIMENT",
                is_fingerprinted=False,
                metadata=None,
                target=target,
                compound_alias=meta["Compound code"],
                structure_origin="EXPERIMENT",
            )

            def add_tag(key, type):
                tagtype, _ = TagType.get_or_create(name=type)
                tag, created = Tag.get_or_create(name=meta[key], type=tagtype)
                if created:
                    tag.clean_and_save()
                pose.tags.add(tag)

            for key, type in POSE_META_TAG_FIELDS.items():
                add_tag(key, type)

            poses.append(pose)

        mrich.success("Loaded", len(poses), "poses for target:", target)

    ### METHODS

    def summary(self) -> None:
        """Print a text summary of this HIPPO"""
        mrich.header(self)
        mrich.var("name", self.name)

        self.db.summary()

        # mrich.var("#compounds", self.num_compounds)
        # mrich.var("#poses", self.num_poses)
        # mrich.var("#targets", self.num_targets)
        # mrich.var("#reactions", self.num_reactions)
        # mrich.var("#tags", self.num_tags)
        # mrich.var("tags", self.tags.unique)

    def get_by_shorthand(self, key) -> "Compound | Pose | Reaction":
        """Get a :class:`.Compound`, :class:`.Pose`, or :class:`.Reaction` by its ID

        :param key: shortname of the object, e.g. C100 for :class:`.Compound` with id=100
        :returns: :class:`.Compound`, :class:`.Pose`, or :class:`.Reaction` object
        """

        if not isinstance(key, str):
            ValueError(f"{key=} must be a str")

        if len(key) > 1:
            ValueError(f"{key=} must be longer than 1")

        prefix = key[0]
        index = key[1:]

        if prefix not in "CPRT":
            raise ValueError(f"Unknown {prefix=}")

        try:
            index = int(index)
        except ValueError:
            mrich.error(f"Cannot convert {index} to integer")
            return None

        match key[0]:
            case "C":
                table = self.compounds
            case "P":
                table = self.poses
            case "T":
                table = self.targets
            case "R":
                table = self.reactions
            case _:
                mrich.error(f"Unsupported {prefix=}")
                return None

        try:
            return table[index]
        except ObjectDoesNotExist:
            mrich.error(f"No object with shorthand: {key=}")

        return None

    ### SINGLE INSERTION

    def register_target(
        self, *, name: str, metadata: str | dict | None = None, debug: bool = False
    ) -> "Target":

        # convert metadata
        from .orm.formatters import dict_formatter

        metadata = dict_formatter(metadata)

        from .protein import Target

        instance, created = Target.get_or_create(name=name, metadata=metadata)

        if created:
            instance.clean_and_save()

            if debug:
                mrich.debug("Created", instance)

        elif debug:
            mrich.debug("Retrieved", instance)

        return instance

    def register_compound(
        self,
        *,
        smiles: str,
        alias: str | None = None,
        metadata: str | dict | None = None,
        debug: bool = False,
    ) -> "Target":

        # convert metadata
        from .orm.formatters import dict_formatter

        metadata = dict_formatter(metadata)

        # flatten molecule
        from .tools import inchikey_from_smiles, sanitise_smiles  # , SanitisationError

        smiles = sanitise_smiles(smiles, verbosity=debug)
        inchikey = inchikey_from_smiles(smiles)

        from .compound import Compound

        # get or create instance
        instance, created = Compound.get_or_create(
            smiles=smiles,
            inchikey=inchikey,
            alias=alias,
            metadata=metadata,
        )

        if created:
            instance.clean_and_save()

            if debug:
                mrich.debug("Created", instance)

        elif debug:
            mrich.debug("Retrieved", instance)

        return instance

    def register_feature(
        self,
        family: str,
        target: "Target | int",
        chain_name: str,
        residue_name: str,
        residue_number: int,
        atom_names: str | dict,
        debug: bool = False,
    ) -> "Feature":

        # convert atom names
        from .orm.formatters import list_formatter

        atom_names = list_formatter(atom_names)

        from .interactions import Feature

        # get related
        from .protein import Target

        if isinstance(target, int):
            target = Target.get(id=target)

        instance, created = Feature.get_or_create(
            family=family,
            target=target,
            chain_name=chain_name,
            residue_name=residue_name,
            residue_number=residue_number,
            atom_names=atom_names,
        )

        if created:
            instance.clean_and_save()

            if debug:
                mrich.debug("Created", instance)

        elif debug:
            mrich.debug("Retrieved", instance)

        return instance

    def register_subsite(
        self,
        name: str,
        target: "Target | int",
        description: str | None = None,
        metadata: str | dict | None = None,
        debug: bool = False,
    ) -> "Subsite":

        # format values
        from .orm.formatters import dict_formatter  # , str_formatter

        metadata = dict_formatter(metadata)
        # description = str_formatter(metadata)

        from .annotation import Subsite

        # get related
        from .protein import Target

        if isinstance(target, int):
            target = Target.get(id=target)

        instance, created = Subsite.get_or_create(
            name=name,
            target=target,
            metadata=metadata,
            description=description,
        )

        if created:
            instance.clean_and_save()

            if debug:
                mrich.debug("Created", instance)

        elif debug:
            mrich.debug("Retrieved", instance)

        return instance

    def register_structure(
        self,
        *,
        protein_file: "str | Path",
        target: "Target | int",
        # is_experimental: bool | None,
        alias: str | None = None,
        pdbblock: str | None = None,
        resolution: float | None = None,
        metadata: str | dict | None = None,
        origin: str = "COMPUTED",
        debug: bool = False,
    ) -> "Structure":

        from .orm.formatters import path_formatter, dict_formatter  # , str_formatter
        from .protein import Target, Structure
        from .files import File, guess_file_format

        # format values
        protein_file = path_formatter(protein_file)
        metadata = dict_formatter(metadata)

        # get related
        if isinstance(target, int):
            target = Target.get(id=target)

        assert origin in Structure.STRUCTURE_ORIGINS

        instance, created = Structure.get_or_create(
            # protein_file=protein_file,
            # is_experimental=is_experimental,
            alias=alias,
            pdbblock=pdbblock,
            resolution=resolution,
            metadata=metadata,
            target=target,
            origin=origin,
        )

        if created:
            instance.clean_and_save()

            if debug:
                mrich.debug("Created", instance)

        elif debug:
            mrich.debug("Retrieved", instance)

        file, file_created = File.get_or_create(
            path=protein_file,
            purpose="INPUT",
            content_type="PROTEIN",
            format_type=guess_file_format(protein_file),
        )

        if created:
            instance.files.add(file)

        return instance

    def register_placement(
        self,
        *,
        # path: "str | Path",
        structure: "Structure | int",
        compound: "Compound | int",
        pose: "Pose | int",
        origin: str = "COMPUTED",
        method: str = "",
        # is_experimental: bool | None,
        # alias: str | None = None,
        # pdbblock: str | None = None,
        # resolution: float | None = None,
        metadata: str | dict | None = None,
        debug: bool = False,
    ) -> "Placement":

        from .orm.formatters import path_formatter, dict_formatter  # , str_formatter
        from .annotation import Placement
        from .protein import Structure
        from .compound import Compound
        from .pose import Pose

        # format values
        # path = path_formatter(path)
        metadata = dict_formatter(metadata)

        # get related
        if isinstance(structure, int):
            structure = Structure.get(id=structure)
        if isinstance(compound, int):
            compound = Compound.get(id=compound)
        if isinstance(pose, int):
            pose = Pose.get(id=pose)

        instance, created = Placement.get_or_create(
            # path=path,
            # is_experimental=is_experimental,
            # alias=alias,
            # pdbblock=pdbblock,
            # resolution=resolution,
            metadata=metadata,
            method=method,
            pose=pose,
            structure=structure,
            compound=compound,
            # target=target,
        )

        if created:
            instance.clean_and_save()

            if debug:
                mrich.debug("Created", instance)

        elif debug:
            mrich.debug("Retrieved", instance)

        return instance

    def register_pose(
        self,
        *,
        smiles: str | None = None,
        mol: "rdkit.Chem.Mol | None" = None,
        complex_file: "str | Path | None" = None,
        protein_file: "str | Path | None" = None,
        mol_file: "str | Path | None" = None,
        alias: str | None = None,
        origin: str = "COMPUTED",
        # is_experimental: bool | None = None,
        is_fingerprinted: bool | None = False,
        metadata: dict | None = None,
        target: "str | Target | int | None" = None,
        target_metadata: dict | None = None,
        compound: "Compound | int | None" = None,
        compound_alias: "str | None" = None,
        compound_metadata: "dict | None" = None,
        placement: "PlacementAttempt | int | None" = None,
        placement_metadata: dict | None = None,
        placement_method: str | None = None,
        structure: "Structure | int | None" = None,
        structure_resolution: float | None = None,
        structure_alias: str | None = None,
        structure_metadata: dict | None = None,
        structure_from_pose: "Pose | int | None" = None,
        structure_origin: str | None = None,
        debug: bool = False,
    ) -> "Pose":
        """Create Pose and optionally related models: PlacementAttempt, Structure"""

        from rdkit import Chem
        from rdkit.Chem.inchi import MolToInchiKey

        from pathlib import Path
        from .tools import sanitise_smiles
        from .orm.formatters import path_formatter, dict_formatter

        from .protein import Target, Structure
        from .compound import Compound
        from .annotation import Placement
        from .pose import Pose
        from .files import File, guess_file_format

        # sanitise / format

        # get related
        if isinstance(target, int):
            target = Target.get(id=target)
        elif isinstance(target, str):
            target = self.register_target(name=target, metadata=target_metadata)

        if isinstance(compound, int):
            compound = Compound.get(id=compound)

        if isinstance(structure, int):
            structure = Structure.get(id=structure)

        if isinstance(placement, int):
            placement = Placement.get(id=placement)

        structure_origin = structure_origin or origin
        assert origin in Pose.POSE_ORIGINS
        assert structure_origin in Pose.POSE_ORIGINS

        protein_file = path_formatter(protein_file)
        complex_file = path_formatter(complex_file)
        mol_file = path_formatter(mol_file)

        # parse molecular files

        if mol:
            if not isinstance(mol, Chem.Mol):
                # convert from binary molecule
                mol = mol.removeprefix(b"MOL\x00")
                mol = Chem.Mol(mol)

            protein_sys = None

        elif mol_file:

            # READ MOL FILE
            # mol =
            raise NotImplementedError

        elif complex_file:

            assert smiles

            # READ COMPLEX FILE
            # mol =
            raise NotImplementedError
            protein_sys = mp.parse(complex_file, verbosity=0).protein_system

        else:

            raise ValueError("Must provide either mol, mol_file, or complex_file")

        if structure:

            # USE STRUCTURE
            protein_file = structure.path
            raise NotImplementedError

        if not protein_sys:

            if complex_file:
                # READ COMPLEX FILE
                protein_sys = mp.parse(complex_file, verbosity=0).protein_system

            elif protein_file:
                # READ PROTEIN FILE
                protein_sys = mp.parse(protein_file, verbosity=0).protein_system

            else:
                raise ValueError(
                    "Must provide either structure, protein_file, or complex_file"
                )

        if not protein_file:
            # look for nearby
            protein_file = complex_file.parent / complex_file.name.replace(
                ".pdb", "_apo-desolv.pdb"
            )
            if not protein_file.exists():
                # write
                protein_sys.write(protein_file, verbosity=0)

        smiles = Chem.MolToSmiles(mol)
        inchikey = MolToInchiKey(mol)

        protein_pdbblock = protein_sys.pdb_block

        # create related (if needed)

        compound = compound or self.register_compound(
            smiles=smiles,
            alias=compound_alias,
            metadata=compound_metadata,
            debug=debug,
        )

        structure = structure or self.register_structure(
            target=target,
            alias=structure_alias,
            resolution=structure_resolution,
            metadata=structure_metadata,
            # is_experimental=is_experimental,
            protein_file=protein_file,
            pdbblock=protein_pdbblock,
            origin=structure_origin,
            debug=debug,
        )

        # placement = placement or self.register_placement(method=placement_method, metadata=placement_metadata)

        # assert related
        assert isinstance(target, Target), (target, type(target))
        assert isinstance(compound, Compound), (compound, type(compound))
        assert isinstance(structure, Structure), (structure, type(structure))

        # raise NotImplementedError

        # format values
        complex_file = path_formatter(complex_file)
        mol_file = path_formatter(mol_file)
        metadata = dict_formatter(metadata)

        pose, created = Pose.get_or_create(
            mol=mol,
            alias=alias,
            smiles=smiles,
            inchikey=inchikey,
            # is_experimental=is_experimental,
            is_fingerprinted=is_fingerprinted,
            metadata=metadata,
            origin=origin,
        )

        if created:
            pose.clean_and_save()

            if debug:
                mrich.debug("Created", pose)

        elif debug:
            mrich.debug("Retrieved", pose)

        placement = placement or self.register_placement(
            compound=compound,
            pose=pose,
            structure=structure,
            metadata=placement_metadata,
            method=placement_method,
            debug=debug,
        )

        if complex_file:
            file, created = File.get_or_create(
                path=complex_file,
                purpose="INPUT",
                content_type="COMPLEX",
                format_type=guess_file_format(complex_file),
            )
            if debug and created:
                mrich.debug("Created:", file)

            pose.files.add(file)

        if protein_file:
            file, created = File.get_or_create(
                path=protein_file,
                purpose="INPUT",
                content_type="PROTEIN",
                format_type=guess_file_format(protein_file),
            )
            if debug and created:
                mrich.debug("Created:", file)

            pose.files.add(file)

        if mol_file:
            file, created = File.get_or_create(
                path=mol_file,
                purpose="INPUT",
                content_type="LIGAND",
                format_type=guess_file_format(mol_file),
            )
            if debug and created:
                mrich.debug("Created:", file)

            pose.files.add(file)

        return pose

    ### BULK INSERTION

    def register_compounds(
        self,
        smiles: list[str],
        alias: str | None = None,
        metadata: str | dict | None = None,
        radical: str = "remove",
        return_objects: bool = True,
        debug: bool = False,
    ) -> "Target":

        from .compound import Compound
        from .tools import inchikey_from_smiles, sanitise_smiles  # , SanitisationError
        from .orm.formatters import dict_formatter

        count_before = len(self.compounds)

        # get or create instance
        objects = []
        inchikeys = []
        for s, a, m in mrich.track(
            zip(smiles, alias, metadata), prefix="sanitising input", total=len(smiles)
        ):

            s = sanitise_smiles(s, verbosity=debug, radical=radical)
            i = inchikey_from_smiles(s)
            m = dict_formatter(m)

            instance = Compound(
                smiles=s,
                inchikey=i,
                alias=a,
                metadata=m,
            )

            objects.append(instance)
            inchikeys.append(i)

        with mrich.loading("Bulk insertion..."):
            created_objs = Compound.bulk_create(
                objects,
                ignore_conflicts=True,
                update_fields=["alias", "metadata"],
                unique_fields=["inchikey", "alias"],
            )

        if return_objects:

            with mrich.loading("Querying ID's..."):
                queried_objs = Compound.filter(inchikey__in=inchikeys)

            with mrich.loading("Mapping inchikeys..."):
                inchikey_map = {obj.inchikey: obj.id for obj in queried_objs}

            with mrich.loading("Assigning ID's..."):
                for obj in created_objs:
                    obj._id = inchikey_map[obj.inchikey]

        num_created = len(self.compounds) - count_before

        if num_created:
            mrich.success("Registered", num_created, "new Compounds")
        else:
            mrich.warning("Registered", num_created, "new Compounds")

        if return_objects:
            return created_objs

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
        pattern = r"[A-Z]\d+"
        if re.match(pattern, key):
            return self.get_by_shorthand(key)

    def __getattr__(self, key: str):
        """Get a :class:`.Compound`, :class:`.Pose`, or :class:`.Reaction` by its ID. See :meth:`.HIPPO.get_by_shorthand`"""
        pattern = r"[A-Z]\d+"
        if re.match(pattern, key):
            return self.get_by_shorthand(key)
