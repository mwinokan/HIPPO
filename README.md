<img src="https://github.com/mwinokan/HIPPO/blob/main/logos/hippo_logo-05.png?raw=true" width="300">

HIPPO
=====

> 🦛 Hit Interaction Profiling for Progression Optimisation

HIPPO is in active development and feedback is appreciated.

Please see the [documentation](https://hippo-docs.winokan.com) to get started


![GitHub Tag](https://img.shields.io/github/v/tag/mwinokan/hippo?include_prereleases&label=PyPI&link=https%3A%2F%2Fpypi.org%2Fproject%2Fhippo-db%2F)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/mwinokan/HIPPO/python-publish.yml?label=publish&link=https%3A%2F%2Fgithub.com%2Fmwinokan%2FHIPPO%2Factions%2Fworkflows%2Fpython-publish.yml)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/mwinokan/HIPPO/black.yml?label=lint&link=https%3A%2F%2Fgithub.com%2Fmwinokan%2FHIPPO%2Factions%2Fworkflows%2Fblack.yml)
[![Documentation Status](https://readthedocs.org/projects/hippo-db/badge/?version=latest)](https://hippo-docs.winokan.com/en/latest/?badge=latest)
![GitHub last commit](https://img.shields.io/github/last-commit/mwinokan/hippo)
![GitHub Issues or Pull Requests](https://img.shields.io/github/issues/mwinokan/hippo)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## Installation

HIPPO is pip-installable, but use of a `conda` environment is recommended for the rdkit and chemicalite dependencies:

```
pip install --upgrade hippo-db
conda install -c conda-forge chemicalite=2024.05.1
```

For more information see the [installation guide](https://hippo-docs.winokan.com/en/latest/#installation)

You can verify the installation:

```
python -m hippo verify
```

Or by running the full suite of tests (see Developer information)

## More Information

<details>

<summary> Information for HIPPO developers </summary>

### Developer installation

To develop on HIPPO please fork this repository and then install locally:

```
git clone https://github.com/YOUR_USER/HIPPO
cd HIPPO
pip install -e .
```

### Releases

HIPPO is automatically released to [PyPI](https://pypi.org/project/hippo-db/) as `hippo-db` via Github [releases](https://github.com/mwinokan/HIPPO/releases) off the `main` branch using the [python-publish](https://github.com/mwinokan/HIPPO/actions/workflows/python-publish.yml) workflow.

### Code style

HIPPO is linted using [black](https://pypi.org/project/black/) and commits are automatically linted using the [black](https://github.com/mwinokan/HIPPO/actions/workflows/black.yml) workflow. The use of [pre-commit](https://pre-commit.com/) is encouraged for local development to automatically run the linting at git commit time:

```
pip install pre-commit
pre-commit install
```

### Documentation

Documentation is automatically built off the [HIPPO/main](https://github.com/mwinokan/HIPPO/tree/main) branch using readthedocs. For local building using sphinx:

```
cd docs
make html
```

To check API reference coverage use [docstr-coverage](https://pypi.org/project/docstr-coverage/)

```
pip install docstr-coverage
docstr-coverage hippo
```

### Tests

Some tests are provided in the tests directory, which can be run with pytest:

```
cd tests
pytest
```

N.B. the numbered tests, e.g. `test_00_cleanup.py` need to run in sequential order to set up the database. Other tests can run in arbitrary order thereafter. The tests will fail if https://fragalysis.diamond.ac.uk can not provide the protein target's data, as specified in tests/config.py.

### Branches

- [HIPPO/main](https://github.com/mwinokan/HIPPO/tree/main): latest stable version
- [HIPPO/dev](https://github.com/mwinokan/HIPPO/tree/dev): @mwinokan's development branch
- [HIPPO/django_lean](https://github.com/mwinokan/HIPPO/tree/django_lean): An experimental branch implementing HIPPO as a Django web-app

</details>
