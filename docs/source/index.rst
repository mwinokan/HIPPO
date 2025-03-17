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

Getting started
===============

HIPPO uses an sqlite database with several inter-connected tables and Python-class representations thereof, the core concepts are explained in :doc:`definitions`. Once familiar you can try :doc:`getting_started`.


.. toctree::
   :maxdepth: 1
   :caption: Documentation Pages

   Home <self>

   Definitions, units, and data types <definitions>

   Getting started <getting_started>

   Adding data <insert_elaborations>

   Interfacing with Syndirella <syndirella>
   
   Running an FFF campaign <fff>
   
   Preparing files for Fragalysis upload <fragalysis>

   Windows installation <windows>

   Random recipe generation <rgen>

   API Reference <api_reference>

.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

.. Structure-based searching <queries>
.. Inserting synthetic pathways <insert_reactions>
.. Quoting with Pycule <pycule_tutorial>