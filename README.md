![hippopotamus_1f99b](https://github.com/mwinokan/HIPPO/assets/36866506/90381bb6-2832-4b12-856f-b8878dfc0dd2)

HIPPO
=====

> 🦛 Hit Interaction Profiling for Procurement Optimisation

A pipeline for optimally selecting building blocks to maximise interaction diversity of reaction products for follow-up screens at XChem.

## Inputs

* Crystallographic data from a fragment screen (i.e. from Fragalysis)
* A set of follow-up compounds with associated catalogue building blocks
* A budget for building block procurement

## Outputs

* [x] Interaction fingerprints for each hit and compound
* [x] Scores for interaction coverage of a given compound set
* [x] Sankey diagram of fragment screen interaction conservation
* [ ] UMAP reduction of the interaction fingerprints into 1D and 2D interaction space
* [ ] Suggested building block sets that optimise interaction coverage and stay within the budget

# Installation

Python dependencies:

* UMAP
* RDKit
* Pandas
* Numpy
* Plotly
* ASE

Install from python source:

* MPyTools
* MolParse
* HIPPO

# Usage

```
import hippo as hp

# create the hippo
pipeline = hp.HIPPO('project_name')

# protein APO structure
pipeline.add_protein_reference(path=protein_reference)

# get the fragment screen data
pipeline.add_hit_metadata(path=metadata_csv)
pipeline.add_hit_directory(path=aligned_path)

# add elaborations/merges
pipeline.add_product_compounds('compounds', metadata_csv, compound_directory, compound_mol_pattern)

# make the interaction fingerprints
pipeline.generate_fingerprints()

# interaction coverage
pipeline.score_interaction_coverage()

# UMAP visualisation (not implemented yet)
pipeline.plot_umap()

# optimise building block selection
pipeline.suggest_building_blocks(budget=10_000)

# store the whole pipeline as a binary
pipeline.write_pickle(pickle_path)

# load an existing pickled pipeline:
pipeline = hp.HIPPO.from_pickle(pickle_path)

```
