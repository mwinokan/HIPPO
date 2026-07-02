.. HIPPO documentation master file, created by
   sphinx-quickstart on Wed Mar 13 16:27:26 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

===================
HIPPO Documentation
===================

*Hit Interaction Profiling for Progression Optimisation* (HIPPO) is a chemical database and python toolkit to expedite fragment-based drug discovery.

Installation
============

Install from PyPI:

::

   $ pip install --upgrade xchem-hippo

For local development, HIPPO uses Docker Compose to provide PostgreSQL with the RDKit cartridge. See the project README for Docker setup instructions.

SQLite is also supported for lightweight local use (no Docker required).

See also :ref:`installation-snippets`.

Getting started
===============

HIPPO uses a database (PostgreSQL with RDKit cartridge, or SQLite for local development) with several inter-connected tables and Python-class representations thereof, the core concepts are explained in :doc:`definitions`. Once familiar you can try :doc:`getting_started`.

.. note::

   HIPPO is built primarily as a Python API to be used in interactive notebooks such as in JupyterLab.

.. toctree::
   :maxdepth: 1
   :caption: Documentation Pages

   Definitions, units, and data types <definitions>

   Getting started <getting_started>

   Adding data <add_data>

   Interfacing with Syndirella and merging algorithms <syndirella>

   API Reference <api_reference>


Core concepts
=============

HIPPO uses a database with several inter-connected tables (see :doc:`db`). In both the database and the python API the following core objects are defined:

Compound
--------

A :class:`.Compound` represents a ligand/small molecule with stereochemistry removed and no atomic coordinates. I.e. it represents the chemical structure. It's name is always an InChiKey. If a compound is an elaboration it can have a :meth:`.Compound.scaffold` property which is another :class:`.Compound`. :class:`.Compound` objects are target-agnostic and can be linked to any number of catalogue entries (:class:`.Quote`) or synthetic pathways (:class:`.Reaction`).

Pose
----

A :class:`.Pose` is a particular conformer of a :class:`.Compound` within a protein environment. A pose will have its own (stereochemical) smiles string, and must have a path to a coordinate file. Poses can have *inspirations* that can be used to trace fragment-derived scaffolds in merges and expansions.

Reaction
--------

A :class:`.Reaction` is a simplified representation of a synthetic pathway to create a product :class:`.Compound`. Reactants (also :class:`.Compound` objects) as well as a reaction type are required.

See :doc:`definitions` for more detail.

.. _installation-snippets:

Installation Snippets
---------------------

If the above fails in your existing software environments, try this:

::

   mamba create --name py312 python=3.12
   mamba activate py312
   pip install xchem-hippo

Additionally, this Dockerfile can be used to create a container with a Jupyter Notebook server:

::

   FROM quay.io/jupyter/minimal-notebook:2025-04-14
   LABEL authors="Max Winokan"

   # Upgrade pip and install packages
   RUN pip install --upgrade pip && pip install xchem-hippo

   # patch rich
   RUN python -c "import mrich; mrich.patch_rich_jupyter_margins()"

   # notebooks
   USER 0
   RUN mkdir "/home/code" && chown ${NB_USER} "/home/code" && \
       sudo apt update && sudo apt install screen -y
   USER ${NB_USER}

   # HIPPO dev branch
   WORKDIR "/home/code"
   RUN git clone https://github.com/mwinokan/HIPPO
   WORKDIR "/home/code/HIPPO"
   RUN git checkout dev && pip install -e . --no-deps

   EXPOSE 8888

   WORKDIR "/home/${NB_USER}"


.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

.. Structure-based searching <queries>
.. Inserting synthetic pathways <insert_reactions>
.. Quoting with Pycule <pycule_tutorial>
