import mrich
from hippo.custom_models import *

from rdkit import Chem
from rdkit.Chem.inchi import MolToInchiKey

from pathlib import Path
from hippo.tools import inchikey_from_smiles, sanitise_smiles
from hippo.orm.formatters import path_formatter, dict_formatter
from hippo.resources import guess_file_format


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
