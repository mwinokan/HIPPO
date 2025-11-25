
======================
Adding data into HIPPO
======================

In theory the :class:`.HIPPO` class makes it simple to insert new material into the database, but there are nuances, especially to improve performance.

In :doc:`getting_started` we saw the :meth:`.HIPPO.add_hits` method which loads all the data from a Fragalysis download or XChemAlign alignment. In this documentation page we'll show how to insert data using other means.

Registering compounds
=====================

The :meth:`.HIPPO.register_compound` method takes as a minimum a smiles string and adds a compound to the database:

::

	c = animal.register_compound(smiles='Clc1c[nH]cn1')

*N.B. to prevent accidents, positional arguments are disabled for many methods in HIPPO.*

In the case that a compound with this structure already exists, the compound is instead retrieved from the database (duplicates will never be created).

Tags, scaffolds, and metadata
-----------------------------

Additional detail can be added to the compound upon registration:

* You may wish to tag the compound with any number of strings.
* You can assign a 'scaffold' compound of which this compound is a superstructure.
* Any dictionary of additional metadata can be stored

::

	c = animal.register_compound(
		smiles='Clc1c[nH]cn1',
		tags=['fragment', 'hit', 'halogen'],
		metadata={'screening_library':'dsiPoised'},
	)

See also the API reference :meth:`.HIPPO.register_compound`.

Registering compounds in bulk
-----------------------------

The :meth:`.HIPPO.register_compounds` method can insert many compounds at once from a list of smiles

::

	smiles = [...]
	values = animal.register_compounds(smiles=smiles)

Part of the HIPPO compound registration process involves sanitising and flattening of smiles and generating InChI-keys. These returned `values` is a list of `(inchikey, new_smiles)` tuples which are in the same order as the input smiles.

Registering a Target and associated Poses
=========================================

The :class:`.Target` class currently only stores a name, but it is a requirement that each :class:`.Pose` is associated with a single :class:`.Target`.

To create (or retrieve) a target:

::

	a71 = animal.register_target(name="A71EV2A")

Then a pose can be registered:

::

	animal.register_pose(
		compound=c,
		target=a71,
		path=/path/to/bound.pdb, # protein ligand complex PDB
	)

Optional arguments include:

* tags
* metadata
* inspirations (other poses that inspired this one)
* reference (in case your path is to a file containing only a ligand .mol, the reference points to another pose who's protein will be used)

See also the API reference :meth:`.HIPPO.register_target` and :meth:`.HIPPO.register_pose`.

Registering poses in bulk
=========================

The :meth:`.Database.register_poses` method can insert many poses at once from a list of dictionaries:

::

    dicts = [
        dict(
            alias=...,           # string can be None
            reference_id=...,    # reference pose id
            inchikey=...,        # pre-computed inchikey
            smiles=...,          # SMILEs
            path=...,            # path to mol-file on disk, used for uniqueness check, can be a fake path
            compound_id=...,     # Compound database ID
            target_id=...,       # Target database ID
            mol=...,             # rdkit.Chem.Mol
            energy_score=...,    # float, can be None
            distance_score=...,  # float, can be None
            metadata=...,        # dictionary, can be empty
        ),
        ...
    ]

    pose_ids = animal.db.register_poses(dicts=dicts)

    poses = animal.poses[pose_ids]

Loading data from an SDF
========================

HIPPO can load compounds, poses, and their related information from a formatted SDF using the :meth:`.HIPPO.load_sdf` method:

::

	animal.load_sdf(
		target=..., # name of target
		path=..., # path to SDF
	)


Refer to the full API reference for more details: :meth:`.HIPPO.load_sdf`

Registering a reaction
======================

To add a reaction:

::

	animal.register_reaction(
		type="Amidation",
		product=c,
		reactants=[a,b],
	)

The above will register an amidation reaction combining compounds **a** and **b** into **c**.

See also the API reference :meth:`.HIPPO.register_reaction`.

