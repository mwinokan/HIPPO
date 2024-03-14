
======================
Adding data into HIPPO
======================

In theory the :class:`.HIPPO` class makes it simple to insert new material into the database, but there are nuances, especially to improve performance.

Registering compounds
=====================

The :meth:`.HIPPO.register_compound` method takes as a minimum a smiles string and adds a compound to the database:

::

	c = animal.register_compound(smiles='Clc1c[nH]cn1')

*N.B. to prevent accidents, positional arguments are disabled for many methods in HIPPO.*

In the case that a compound with this structure already exists, the compound is instead retrieved from the database (duplicates will never be created).

Tags, bases, and metadata
-------------------------

Additional detail can be added to the compound upon registration:

* You may wish to tag the compound with any number of strings.
* You can assign a 'base' compound of which this compound is a superstructure.
* Any dictionary of additional metadata can be stored

::

	c = animal.register_compound(
		smiles='Clc1c[nH]cn1',
		tags=['fragment', 'hit', 'halogen'],
		metadata={'screening_library':'dsiPoised'},
	)

See also the API reference :meth:`.HIPPO.register_compound`.

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

Improving Performance
=====================

Querying a large database can become expensive, and permance gains can also be realised by commiting database transactions in bulk. Additionally, initialising python object instances like Compound may not be necessary.

In a large loop of registering compounds, reactions, and poses the above issues can be mitigated:

::

	for i,row in df.iterrows():

		# register the reactants

		reactant1 = animal.register_compound(
			smiles=row.r1_smiles, 
			return_compound=False, # this just returns the compound ID
			commit=False, # this does not commit the change to the database
		)

		reactant2 = animal.register_compound(
			smiles=row.r2_smiles, 
			return_compound=False, # this just returns the compound ID
			commit=False, # this does not commit the change to the database
		)

		# register the product

		product = animal.register_compound(
			smiles=row.smiles, 
			return_compound=True, 
			metadata=comp_metadata, 
			commit=False
		)

		# register the reaction

		reaction = animal.register_reaction(
			type=row.reaction, 
			product=product, 
			reactants=[reactant1, reactant2], 
			commit=False,
		)

		# register the pose
		
		pose = animal.register_pose(
			compound=product, 
			target='A71EV2A', 
			path=pose_path, 
			metadata=pose_metadata, 
			inspirations=[inspiration1, inspiration2], 
			commit=False, 
			return_pose=False, 
			overwrite_metadata=True, # don't bother checking existing metadata
		)

		# add tags and bases:

		animal.db.update(
			table='compound', 
			id=product_id, 
			key='compound_base', 
			value=base_id, 
			commit=False
		)

		animal.db.insert_tag(
			name='elab', 
			compound=product_id, 
			commit=False,
		)

		animal.db.commit()
