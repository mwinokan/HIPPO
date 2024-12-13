import mrich


def pose_path_mapper(kwargs) -> str:

    # mrich.print(kwargs)

    path = kwargs["path"]

    del kwargs["path"]

    if path.endswith(".pdb"):
        kwargs["complex_file"] = path

    # mrich.print(kwargs)

    return kwargs


def pose_reference_mapper(kwargs) -> str:

    structure_from_pose = kwargs["structure_from_pose"]

    if not structure_from_pose:
        del kwargs["structure_from_pose"]

    return kwargs


LEGACY_MAP = [
    # targets
    {
        "table": "target",
        "model": "Target",
        "method": "register_target",
        "fields": {
            "target_name": "name",
            "target_metadata": "metadata",
        },
    },
    # # compounds
    # {
    # 	"table": "compound",
    # 	"bulk":True,
    # 	"model": "Compound",
    # 	"method": "register_compounds",
    # 	"fields": {
    # 		"compound_alias":"alias",
    # 		"compound_smiles":"smiles",
    # 		"compound_metadata":"metadata",
    # 	},
    # },
    # poses
    {
        "table": "pose",
        "model": "Pose",
        "method": "register_pose",
        "fields": {
            # "pose_inchikey": "",
            "pose_alias": "alias",
            # "pose_smiles": "",
            "pose_reference": "structure_from_pose",
            "pose_path": "path",
            # "pose_compound": "compound",
            "pose_target": "target",
            "pose_mol": "mol",
            "pose_fingerprint": "is_fingerprinted",
            # "pose_energy_score": "",
            # "pose_distance_score": "",
            "pose_metadata ": "metadata",
        },
        "custom_fields": {
            "pose_path": pose_path_mapper,
            "pose_reference": pose_reference_mapper,
        },
        "related_fields": {
            "structure_from_pose": "pose",
        },
        "warnings": {
            "Not migrating scores",
        },
        "debug": False,
    },
    # features
    {
        "table": "feature",
        "model": "Feature",
        "method": "register_feature",
        "fields": {
            "feature_family": "family",
            "feature_target": "target",
            "feature_chain_name": "chain_name",
            "feature_residue_name": "residue_name",
            "feature_residue_number": "residue_number",
            "feature_atom_names": "atom_names",
        },
        "related_fields": {"target": "target"},
    },
    # subsites
    {
        "table": "subsite",
        "model": "Subsite",
        "method": "register_subsite",
        "fields": {
            "subsite_target": "target",
            "subsite_name": "name",
            "subsite_metadata": "metadata",
        },
        "related_fields": {"subsite_target": "target"},
    },
]
