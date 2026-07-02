
==========================
Getting started with HIPPO
==========================

To create a HIPPO animal, use the :func:`~hippo.bootstrap.load_hippo` factory function (exposed as ``HIPPO``):

::

	from hippo import HIPPO

	animal = HIPPO(
		target_name="A71EV2A",
		target_access_string="lb18145-1",
		username="your_fed_id",
	)

This configures Django, connects to the database, and returns an animal object bound to the specified target. Connection parameters are read from environment variables (``DB_NAME``, ``DB_USER``, ``DB_PASSWORD``, ``DB_HOST``, ``POSTGRES_PORT``) or can be passed explicitly via the ``db`` parameter.

For local SQLite usage (no Docker needed):

::

	animal = HIPPO(
		target_name="A71EV2A",
		target_access_string="lb18145-1",
		username="your_fed_id",
		db="path/to/db.sqlite",
	)

.. note::

	For PostgreSQL, ensure Docker Compose is running (``docker compose up database``) and your ``.env`` file is configured. See the project README for setup instructions.


Loading crystallographic hits from Fragalysis
=============================================

1. Go `Fragalysis <https://fragalysis.diamond.ac.uk>`_ and download a dataset. Alternatively, you can use the `Fragalysis Python API <https://fragalysis.readthedocs.io/en/latest/py_api.html>`_.

2. Load the crystallographic data

::

	animal.add_hits(
		metadata_csv='/path/to/metadata.csv',
		aligned_directory='/path/to/aligned_files',
	)

The animal is already bound to a target at initialization, so ``target_name`` is not needed here.

.. attention::

	N.B. all poses loaded into a HIPPO database only have an absolute path stored to the original file - they are not copied! It is your responsibility to ensure that their original files remain accessible.


Registering methods
===================

Before loading posed virtual hits, register the computational methods used to generate them:

::

	animal.register_pose_method(name="xray", version="1.0.0", description="Crystallographic poses")
	animal.register_pose_method(name="fragmenstein", version="1.0.0", description="Fragmenstein placement")
	animal.register_scoring_method(name="gnina_cnn_vs", version="1.3.2", description="GNINA CNN VS score")
	animal.register_enumeration_method(name="fragmenstein", version="1.0.0", description="Fragmenstein merges")


Navigating compounds and poses
==============================

The below sections explain how to work with :class:`.Compound` objects and sets thereof. For further details on the concept of :class:`.Compound` objects and others see :doc:`definitions`.

Getting compounds/poses
-----------------------

Compounds can be accessed via the compounds property which returns a :class:`.CompoundSet`:

::

	all_compounds = animal.compounds

CompoundSets can be indexed using a database id (positive integer) or sliced:

::

	c = animal.compounds[1]
	subset = animal.compounds[20:30]

You can filter compounds by tag:

::

	hits = animal.compounds.get_by_tag('hits')

.. See also the :doc:`tools for structure-based searching<queries>`

Equivalent methods exist for animal.poses (returns a :class:`.PoseSet`), animal.reactions (returns a :class:`.ReactionSet`), and animal.interactions (returns a :class:`.InteractionSet`). See also the :doc:`api_reference` pages.

Inspecting a compound and its poses
-----------------------------------

Once you have a compound you can access database properties using its properties:

::

	c = animal.compounds[1]

	c.id # Database ID (int)
	c.name # alias or InChiKey
	c.smiles # (flattened) smiles
	c.mol # rdkit.Chem.Mol
	c.tags # assigned tags
	c.metadata # metadata dictionary

	c.draw() # draw the molecule (and its scaffold)

You can access a compounds poses, which have similar functionality

::

	poses = c.poses

	p = poses[0]

	p.id # Database ID (int)
	p.pose_alias # pose name
	p.pose_smiles # (stereo) smiles
	p.mol # rdkit.Chem.Mol
	p.pose_metadata # metadata dictionary

	p.draw() # draw the molecule pose (3d)

See also the API reference for :doc:`compounds <compounds>` and :doc:`poses <poses>`.

Interaction fingerprinting
==========================

Interactions are fingerprinted as follows:

::

	import mrich

	for pose in mrich.track(animal.poses):
		pose.calculate_interactions()

.. note::

	``mrich.track`` just gives you a nice progress bar.

Graphing
========

Interaction fingerprints can be visualised with a punchcard:

::

	animal.plot_interaction_punchcard(poses=animal.poses(tag='hits'), subtitle='hits', group='pose_name')

See also :func:`.plotting.plot_interaction_punchcard`.
