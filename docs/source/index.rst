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

Hit Interaction Profiling for Procurement Optimisation is a chemical database and python toolkit to make informed sampling decisions for efficient SAR exploration.

HIPPO is still in alpha-development, but it can be obtained from PyPI:

::
   
   pip install --upgrade hippo-db

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


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   Getting started <getting_started>

   Structure-based searching <queries>

   Inserting virtual hits <insert_elaborations>
   
   Inserting synthetic pathways <insert_reactions>

   Quoting with Pycule <pycule_tutorial>

   API Reference <api_reference>


.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
