
==========================
Loading Syndirella outputs
==========================

Syndirella has been developed to produce a HIPPO-friendly pickled SDF output following the syntax:

::

	{inchikey}_{reaction_uuid}_to_hippo.pkl.gz

Loading base compounds/compound scaffolds
=========================================

If a single pickled DataFrame is available of all base compounds, then the method :meth:`.HIPPO.add_syndirella_bases` can be used. However, as this is not always the case this page documents using the :meth:`.HIPPO.add_syndirella_elabs` method with the ``base_only=True`` keyword argument:

::

	animal.add_syndirella_elabs(path, base_only=True)

Using the inspiration_map argument
----------------------------------

Often the inspiration fragments are referenced using the observation long-code, in order to map to a :prop:`.Pose.id` in a performant manner an ``inspiration_map`` can be used. For example it could be a dictionary:

::

	inspiration_map = {
		"CHIKV_MacB-x1203_A_501_CHIKV_MacB-x0300+A+401+1": 10,
		...
	}

If there is a reliable syntax, a function can be passed instead:

::

	def inspiration_map(longcode):
	    crystal_dataset = longcode.removeprefix('CHIKV_MacB-')[:5]
	    return animal.poses[f'c{crystal_dataset}a'].id

