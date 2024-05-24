## Setup

1. Fragalysis download of crystallographic data
1. Create blank hippo DB
1. Load aligned hits into DB
1. Register an extra pose for Ax0310a (Matteo's relaxed) (store pose_ID somewhere global)

## Base compounds

1.  Get inspiration 1 `inspiration1 = animal.poses[...]`
1.  Get inspiration 2 `inspiration2 = animal.poses[...]`
1.  Register the base compound
    1. `tags = ['base']`
    2. `metadata = dict(warren_name=..., )`
    1. `base = animal.register_compound(smiles=smiles, tags=tags, metadata=metadata)`

1. Place the compound with Fragmenstein
2. Register a pose
   1. `tags = ['base', 'flat_placement']`
   1. `metadata = dict(ddG=..., mRMSD=...)`
   1. `pose = animal.register_pose(smiles=smiles, compound=base, tags=tags, metadata=metadata, inspirations=[inspiration1, inspiration2])`

## Elaborations & placement pipeline (parralel jobs)

1. same as before but stop at num_atom_diff == 15

1. Write a new CSV output with columns:
   * compound_ID (from hippo DB)
   * paths to minimised .mol (only if successful?)
   * reactant smiles
   * reaction names
   * ddG
   * RMSD
   * acceptable?
   * errors?

## Load in poses (single thread / job)

1. Parse through the new placement CSV's
1. Register poses if ddG, RMSD, acceptable, errors columns are all OK
1. For each registered pose add the following detail:
   * tags = ['base'] or ['elab']
   * metadata = dict(ddG=..., RMSD=...)
