
=====================
Definitions and units
=====================

Definitions
===========

HIPPO uses an sqlite database with several inter-connected tables (see :doc:`db`). In both the database and the python API the following core objects are defined:

Target
------

A :class:`.Target` represents a certain protein / XCA-alignment as uploaded to Fragalysis.

Compound
--------

A :class:`.Compound` represents a ligand/small molecule with stereochemistry removed and no atomic coordinates. I.e. it represents the chemical structure. It's default name is always an InChiKey. :class:`.Compound` objects can have an :attr:`.Compound.alias` which is a custom name which will supercede the InChiKey when representing the compound. :class:`.Compound` objects also have a shorthand prefixed with ``C``, for example: ``C1`` which refers to the compound with database id 7273.

:: 

	c1 = animal.register_compound(smiles="OCc1ccc2c(c1)CCO2")
	print(c1)
	c1.draw()

.. image:: ../images/compound_output.png
  :width: 381
  :alt: compound_output

.. seealso::
	- :attr:`.HIPPO.compounds`
	- :class:`.Compound`
	- :class:`.CompoundSet`
	- :class:`.CompoundTable`

Bases / Elaborations
~~~~~~~~~~~~~~~~~~~~

Scaffold / superstructure relationships can also be encoded for :class:`.Compound` objects. Namely, the :attr:`.Compound.bases` property can be used to access other :class:`.Compound` objects that have been labelled as scaffolds/bases/substructures, and :attr:`.Compound.elabs` is used to access the inverse relationship.

:: 

	c2 = animal.register_compound(smiles="OCc1ccc2c(c1F)CCO2")
	c2.add_base(base=c1)
	print(c2)
	c2.draw()

.. image:: ../images/compound_elab.png
  :width: 403
  :alt: compound_elab

.. seealso::
	- :attr:`.HIPPO.bases`
	- :attr:`.HIPPO.elabs`
	- :attr:`.Compound.bases`
	- :attr:`.Compound.elabs`
	- :attr:`.CompoundSet.bases`
	- :attr:`.CompoundSet.elabs`
	- :attr:`.CompoundTable.bases`
	- :attr:`.CompoundTable.elabs`

Ingredient
----------

An :class:`.Ingredient` is defined as a specific quantity (in ``mg``) of a :class:`.Compound` and is used when defining quotes and recipes.

.. seealso::
	- :class:`.Ingredient`
	- :class:`.IngredientSet`

Pose
----

A :class:`.Pose` is a particular conformer of a :class:`.Compound` within a protein environment. A pose will have its own (stereochemical) smiles string, and must have a path to a coordinate file. This file can either be a ``.mol`` molecule file or a ``.pdb`` file of the protein-ligand complex.

:: 

	p1 = c1.poses[0]
	print(p1)
	p1.draw()

.. image:: ../images/pose_output.png
  :width: 364
  :alt: pose_output

Reference
~~~~~~~~~

When a pose has been defined from a ``.mol`` file without a protein conformation, a :attr:`.Pose.reference` can be set to use the protein conformation from a different pose.

Inspirations
~~~~~~~~~~~~

Poses can have :attr:`.Pose.inspirations` that can be used to link to other poses that have been referenced in the design of this pose, for example it can be used to link to experimental fragment hits referenced in the fragment-growing/merging compound design.

.. seealso::
	- :attr:`.HIPPO.poses`
	- :class:`.Pose`
	- :class:`.PoseSet`
	- :class:`.PoseTable`

Tag
---

:class:`.Compound` and :class:`.Pose` objects can be tagged with arbitrary :class:`.Tag` strings to categorise them.

.. seealso::
	- :attr:`.HIPPO.tags`
	- :attr:`.Compound.tags`
	- :attr:`.Pose.tags`
	- :class:`.TagSet`
	- :class:`.TagTable`

Quote
-----

Procurement and catalogue/inventory availability information for :class:`.Compound` entries can be added to the database and interfaced with :class:`.Quote` objects.

.. seealso::
	- :class:`.Quote`
	- :attr:`.Ingredient.quote`
	- :meth:`.Compound.get_quotes`
	- :meth:`.Ingredient.get_quotes`

Reaction
--------

A :class:`.Reaction` is a simplified representation of a chemical reaction from the :attr:`.Reaction.reactants` (:class:`.CompoundSet`) to a single :attr:`.Reaction.product` (:class:`.Compound`).

.. seealso::
	- :attr:`.HIPPO.reactions`
	- :class:`.Reaction`
	- :class:`.ReactionSet`
	- :class:`.ReactionTable`

Recipe
------

A :class:`.Recipe` describes a synthetic pathway, potentially containing multiple :class:`.Reaction` steps to any number of :class:`.Compound` products (:attr:`Recipe.products`). Recipes are not stored in the database but can be serialized into ``JSON``.

.. seealso::
	- :class:`.Recipe`

	- :class:`.RandomRecipeGenerator`

Route
-----

A :class:`.Route` is a special case of the :class:`.Recipe` mechanism, with the distinction that it encodes the information needed to synthesise a single product :class:`.Compound`. Routes can be stored and retrieved from the database.

.. seealso::
	- :class:`.Route`
	- :class:`.RouteSet`

Subsite
-------

:class:`.Subsite` records are an additional annotation that can be applied to :class:`.Pose` entries, these should be used to indicate which subsites a pose occupies on a protein target.

.. seealso::
	- :attr:`.HIPPO.Subsites`
	- :attr:`.Pose.Subsites`
	- :class:`.Subsite`

Feature
-------

A :class:`.Feature` is a pharmacophoric feature on a given protein :class:`.Target`

.. seealso::
	- :attr:`.Target.features`
	- :class:`.Feature`

Interaction
-----------

The :class:`.Interaction` class can be used to store protein-ligand interactions between pharmacophores on the ligand and :class:`.Feature` records.

.. seealso::
	- :attr:`.HIPPO.interactions`
	- :attr:`.Pose.interactions`
	- :meth:`.Pose.calculate_interactions`
	- :meth:`.Pose.calculate_prolif_interactions`
	- :class:`.Interaction`
	- :class:`.InteractionSet`
	- :class:`.InteractionTable`

Units
=====

- lead time: days
- compound quantities/amounts: mg
- purity: fraction [0,1]
- product_yield: fraction [0,1]
