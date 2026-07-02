
======================
Adding data into HIPPO
======================

In theory the :class:`.HIPPO` class makes it simple to insert new material into the database, but there are nuances, especially to improve performance.

In :doc:`getting_started` we saw the :meth:`.HIPPO.add_hits` method which loads all the data from a Fragalysis download or XChemAlign alignment. In this documentation page we'll show how to insert data using other means.

Registering compounds
=====================

.. note::

	In the current Django-based HIPPO, compounds are created implicitly when loading data via :meth:`.HIPPO.add_hits` or :meth:`.HIPPO.load_sdf`. There is no standalone ``register_compound`` method on the HIPPO class.

Tags, scaffolds, and metadata
-----------------------------

Compounds can be tagged and annotated with metadata at load time. When loading from SDFs, use the ``compound_tags`` parameter of :meth:`.HIPPO.load_sdf`. When loading hits, tags can be passed to :meth:`.HIPPO.add_hits`.

Registering compounds in bulk
-----------------------------

Compounds are registered in bulk implicitly when loading SDFs via :meth:`.HIPPO.load_sdf`. Deduplication is handled automatically based on InChI-Keys.

Registering Poses
=================

.. note::

	In the current Django-based HIPPO, the target is set at initialization time (when calling ``HIPPO(target_name=...)``). There is no need to explicitly register targets.

Poses are typically loaded via :meth:`.HIPPO.add_hits` (for crystallographic data) or :meth:`.HIPPO.load_sdf` (for virtual hits from SDFs). See :doc:`getting_started` for examples.

Loading poses in bulk
=====================

Poses are loaded in bulk via :meth:`.HIPPO.load_sdf` which handles compound creation, pose registration, and deduplication in one step.

Loading data from an SDF
========================

HIPPO can load compounds, poses, and their related information from a formatted SDF using the :meth:`.HIPPO.load_sdf` method:

::

	animal.load_sdf(
		path=..., # path to SDF
	)


Refer to the full API reference for more details: :meth:`.HIPPO.load_sdf`

Registering a reaction
======================

.. note::

	In the current Django-based HIPPO, reactions are typically created via synthesis route loading (e.g. :meth:`.HIPPO.add_syndirella_routes`). There is no standalone ``register_reaction`` method on the HIPPO class.
