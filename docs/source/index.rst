.. HIPPO documentation master file, created by
   sphinx-quickstart on Wed Mar 13 16:27:26 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

===================
HIPPO Documentation
===================

.. image:: ../../logos/hippo_logo-05.png
  :width: 400
  :alt: HIPPO logo

Hit Interaction Profiling for Procurement Optimisation is a chemical database and python toolkit to make informed sampling decisions for effective SAR exploration.

N.B. HIPPO and this documentation are still in alpha-development.

Installation
============

On Mac OS and Linux it is recommended to install from PyPI using Conda/Miniconda. 

Chemicalite is not supported on Windows, but there is a workaround described in the :doc:`windows`.

The `hippo` python module can be obtained from PyPI:

::
   
   $ pip install --upgrade hippo-db

You will also need `chemicalite` which is an extension to SQLite for cheminformatics:

::

   $ conda install -c conda-forge chemicalite=2022.04.1

N.B. Compatibility between rdkit and chemicalite versions is quite strict, and database files created with a certain version pair may not be interoperable with others. 

Core concepts
=============

HIPPO uses an sqlite database with several inter-connected tables (see :doc:`db`). In both the database and the python API the following core objects are defined:

Compound
--------

A :class:`.Compound` represents a ligand/small molecule with stereochemistry removed and no atomic coordinates. I.e. it represents the chemical structure. It's name is always an InChiKey. If a compound is an elaboration it can have a :meth:`.Compound.base` property which is another :class:`.Compound`. :class:`.Compound` objects are target-agnostic and can be linked to any number of catalogue entries (:class:`.Quote`) or synthetic pathways (:class:`.Reaction`). 

Pose
----

A :class:`.Pose` is a particular conformer of a :class:`.Compound` within a protein environment. A pose will have its own (stereochemical) smiles string, and must have a path to a coordinate file. Poses can have *inspirations* that can be used to trace fragment-derived scaffolds in merges and expansions.

Reaction
--------

A :class:`.Reaction` is a simplified representation of a synthetic pathway to create a product :class:`.Compound`. Reactants (also :class:`.Compound` objects) as well as a reaction type are required.

See :doc:`definitions` for more detail.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   Definitions, units, and data types <definitions>

   Getting started <getting_started>

   Structure-based searching <queries>

   Adding data <insert_elaborations>

   Loading Syndirella outputs <load_syndirella>
   
   Inserting synthetic pathways <insert_reactions>

   Preparing files for Fragalysis upload <fragalysis>

   Quoting with Pycule <pycule_tutorial>

   Windows installation <windows>

   Random recipe generation <rgen>

   API Reference <api_reference>

.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
