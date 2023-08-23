![hippopotamus_1f99b](https://github.com/mwinokan/HIPPO/assets/36866506/90381bb6-2832-4b12-856f-b8878dfc0dd2)

HIPPO
=====

> 🦛 Hit Interaction Profiling for Procurement Optimisation

A pipeline for optimally selecting building blocks to maximise interaction diversity of reaction products for follow-up screens at XChem.

## Inputs

* [ ] Crystallographic data from a fragment screen (i.e. from Fragalysis)
* [ ] A set of follow-up compounds with associated catalogue building blocks
* [ ] A budget for building block procurement

## Outputs

* [ ] Interaction fingerprints for each hit and compound
* [ ] Scores for interaction coverage of a given compound set
* [ ] Sankey diagram of fragment screen interaction conservation
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
pipeline = hp.HIPPO()

crystal_path = '...'
compound_paths = [
  ...
]

# assign paths
pipeline.get_crystallographic_directory(crystals=crystal_path)

# load in the compound sets
for path in compound_paths():
  pipeline.add_compound_set(path=path)

# interaction coverage
pipeline.generate_fingerprints()
pipeline.score_interaction_coverage()

# UMAP visualisation
pipeline.plot_umap()

# optimise building block selection
pipeline.suggest_building_blocks(budget=10_000)

```
