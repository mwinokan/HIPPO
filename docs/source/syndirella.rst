
===========================
Interfacing with Syndirella
===========================

`Syndirella <http://github.com/kate-fie/syndirella>`_ is a python package used to enumerate synthetically accessible chemical space around **scaffold** compounds by first performing retrosynthesis calculations, superstructure searches of reactants, cartesian multiplication of reactant superstructure sets, and placement with `Fragmenstein <http://github.com/matteoferla/fragmenstein>`_. Syndirella is used extensively for fragment progression at XChem, and HIPPO has been co-developed to interface with it.

Generating Syndirella inputs
============================

Syndirella requires a CSV input with details for each scaffold. Since scaffolds have to be placed relative to inspiration structures, and into a protein template. The HIPPO :class:`.Pose` object is most suitable as it encodes both :attr:`.Pose.inspirations` and :attr:`.Pose.reference` properties.

Syndirella can be run as a single process which sequentially processes scaffolds, or distributed across SLURM jobs.

Within the specified directory, :meth:`.PoseSet.to_syndirella` will generate:

- a ``templates`` subdirectory with apo-structures for each reference protein
- an input CSV compatible with syndirella with a row for each scaffold
- an SDF containing all inspiration ligands

Processing all scaffolds in a single syndirella run
---------------------------------------------------

.. code-block:: python

    scaffolds = animal.poses[...]
    scaffolds.to_syndirella("elabs")

Creating separate inputs for each scaffold
------------------------------------------

.. code-block:: python

    scaffolds = animal.poses[...]
    for i,pose in enumerate(scaffolds):
        mrich.h3(f"{i}/{len(scaffolds)} {pose.id}")
        pset = animal.poses[pose.id,]
        pset.to_syndirella(f"elabs/P{pose.id}")

Running Syndirella
==================

Syndirella will need to be run on the files generated in the previous section. An example script which will need to be run for each set of inputs specified as a command-line argument:

.. code-block:: bash

    #!/bin/bash

    set -e

    KEY=$1

    echo --input $(pwd)/$KEY"_syndirella_input.csv"
    echo --hits_path $(pwd)/$KEY"_syndirella_inspiration_hits.sdf"
    echo --output $(pwd)/$KEY"_elabs"
    echo --metadata METADATA_CSV_PATH
    echo --templates "$(pwd)/templates"

    /opt/xchem-fragalysis-2/maxwin/conda/bin/syndirella \
        --input $(pwd)/$KEY"_syndirella_input.csv" \
        --hits_path $(pwd)/$KEY"_syndirella_inspiration_hits.sdf" \
        --output $(pwd)/$KEY"_elabs" \
        --templates "$(pwd)/templates" \
        --metadata METADATA_CSV_PATH \
        --no_scaffold_place

The ``METADATA_CSV_PATH`` placeholder will need to be replaced.

Loading Syndirella outputs
==========================

Syndirella has been developed to produce a HIPPO-friendly output in the syntax ``{inchikey}_{reaction_uuid}_to_hippo.pkl.gz``. This can be parsed using the :meth:`.HIPPO.add_syndirella_elabs` method. In the case where many such files are loaded it is recommended to use a shell script run via SLURM, such as the one below:

.. code-block:: python
    :linenos:
    :emphasize-lines: 17
    :caption: Example script for loading syndirella elaborations and their routes

    from pathlib import Path
    import hippo
    import mrich

    mrich.var("hippo", hippo.__file__)

    animal = hippo.HIPPO(PROJECT_NAME, DATABASE_PATH)
    animal.db.backup()

    output_root = Path("../syndirella/elabs/")

    files = list(output_root.glob("*/*-*-?/*to_hippo*"))

    for i,file in enumerate(files):
        mrich.h2(f"{i+1}/{len(files)}")
        try:
            animal.add_syndirella_elabs(file)
        except Exception as e:
            mrich.error(file)
            mrich.error(e)
            continue
        
    animal.db.close()

The above script can be submitted to the DLS / IRIS cluster as follows:

.. code-block:: bash

    sbatch --job-name load_elabs --exclusive --no-requeue /opt/xchem-fragalysis-2/maxwin/slurm/run_python.sh 4_load_elabs.py
