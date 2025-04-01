import mrich
from pathlib import Path
import molparse as mp

from rdkit import Chem
from rdkit.Chem.inchi import MolToInchiKey

from hippo.custom_models import *
from hippo.tools import inchikey_from_smiles, sanitise_smiles
from hippo.orm.formatters import path_formatter, dict_formatter
from hippo.resources import guess_file_format


def register_target(
    *,
    name: str,
    metadata: str | dict | None = None,
    debug: bool = False,
) -> "Target":

    # convert metadata
    metadata = dict_formatter(metadata)

    instance, created = Target.get_or_create(name=name, metadata=metadata)

    if created:
        instance.clean_and_save()

        if debug:
            mrich.debug("Created", instance)

    elif debug:
        mrich.debug("Retrieved", instance)

    return instance


def register_compound(
    *,
    smiles: str,
    alias: str | None = None,
    metadata: str | dict | None = None,
    debug: bool = False,
) -> "Target":

    # convert metadata
    metadata = dict_formatter(metadata)

    # flatten molecule
    smiles = sanitise_smiles(smiles, verbosity=debug)
    inchikey = inchikey_from_smiles(smiles)

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


def register_pose(
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

    # sanitise / format

    # get related
    if isinstance(target, int):
        target = Target.get(id=target)
    elif isinstance(target, str):
        target = register_target(name=target, metadata=target_metadata)

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
            # raise ValueError(
            # "Must provide either structure, protein_file, or complex_file"
            # )
            mrich.warning("No structure, protein_file, or complex_file provided")

    if not protein_file and complex_file:

        # look for nearby
        protein_file = complex_file.parent / complex_file.name.replace(
            ".pdb", "_apo-desolv.pdb"
        )
        if not protein_file.exists():
            # write
            protein_sys.write(protein_file, verbosity=0)

    smiles = Chem.MolToSmiles(mol)
    inchikey = MolToInchiKey(mol)

    if protein_sys:
        protein_pdbblock = protein_sys.pdb_block
    else:
        protein_pdbblock = None

    # create related (if needed)

    compound = compound or register_compound(
        smiles=smiles,
        alias=compound_alias,
        metadata=compound_metadata,
        debug=debug,
    )

    if protein_file and protein_pdbblock:

        structure = structure or register_structure(
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

    else:

        mrich.warning("Skipping Structure registration")

    # placement = placement or register_placement(method=placement_method, metadata=placement_metadata)

    # assert related
    assert isinstance(compound, Compound), (compound, type(compound))
    assert target is None or isinstance(target, Target), (target, type(target))
    assert structure is None or isinstance(structure, Structure), (
        structure,
        type(structure),
    )

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

    if structure:
        placement = placement or register_placement(
            compound=compound,
            pose=pose,
            structure=structure,
            metadata=placement_metadata,
            method=placement_method,
            debug=debug,
        )
    else:
        placement = None

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


def register_structure(
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

    # format values
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


def register_fragalysis_target(
    target_name: str,
    target_access_string: str | None = None,
    download_directory: "str | Path | None" = None,
    metadata_csv: "str | Path | None" = None,
    aligned_dir: "str | Path | None" = None,
    stack: str = "production",
    token: str | None = None,
):

    import pandas as pd
    from rdkit.Chem import PandasTools
    from .fragalysis import download_target, POSE_META_TAG_FIELDS

    if target_access_string:

        from django.conf import settings

        assert download_directory, "download_directory required"

        download_directory = Path(download_directory)

        if not download_directory.exists():
            download_directory.mkdir()

        target_dir = download_target(
            target_name=target_name,
            target_access_string=target_access_string,
            destination=download_directory,
            stack=stack,
            token=token,
        )

        if not target_dir:
            raise ValueError(
                "No target, is Fragalysis available and are permissions correct?"
            )

        metadata_csv = target_dir / "metadata.csv"
        aligned_dir = target_dir / "aligned_files"

    else:

        assert (
            metadata_csv
        ), "metadata.csv must be provided if no target_access_string is given"
        assert (
            aligned_dir
        ), "aligned_dir must be provided if no target_access_string is given"
        metadata_csv = Path(metadata_csv)
        aligned_dir = Path(aligned_dir)

    assert metadata_csv.exists()
    assert aligned_dir.exists()

    meta_df = pd.read_csv(metadata_csv)

    target = register_target(name=target_name)

    shortcodes = set(meta_df["Code"].values)

    poses = []

    subdirs = list(aligned_dir.iterdir())

    for subdir in mrich.track(subdirs):

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

        pose = register_pose(
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
