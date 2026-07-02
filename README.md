<img src="https://github.com/xchem/HIPPO/blob/main/logos/hippo_logo-05.png?raw=true" width="300">

# XChem HIPPO

> HIPPO: 🦛 Hit Interaction Profiling for Progression Optimisation

HIPPO is in active development and feedback is appreciated.

Please see the [documentation](https://hippo-docs.winokan.com) to get started

![GitHub Tag](https://img.shields.io/github/v/tag/xchem/hippo?include_prereleases&label=PyPI&link=https%3A%2F%2Fpypi.org%2Fproject%2Fxchem-hippo%2F)
![Release](https://img.shields.io/github/actions/workflow/status/xchem/HIPPO/release.yaml?label=publish&link=https%3A%2F%2Fgithub.com%2Fxchem%2FHIPPO%2Factions%2Fworkflows%2Frelease.yaml)
![Lint](https://img.shields.io/github/actions/workflow/status/xchem/HIPPO/lint.yaml?label=lint&link=https%3A%2F%2Fgithub.com%2Fxchem%2FHIPPO%2Factions%2Fworkflows%lint.yaml)
![Test](https://img.shields.io/github/actions/workflow/status/xchem/HIPPO/test.yaml?label=test&link=https%3A%2F%2Fgithub.com%2Fxchem%2FHIPPO%2Factions%2Fworkflows%test.yaml)

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## Installation

HIPPO is pip-installable:

```bash
pip install --upgrade xchem-hippo
```

For local development with PostgreSQL + RDKit cartridge, use Docker Compose (see below).

For more information see the [installation guide](https://hippo-docs.winokan.com/en/latest/#installation)

## More Information

<details>

<summary>Repository structure</summary>

### Branches

- [HIPPO/main](https://github.com/xchem/HIPPO/tree/main): latest stable version (Django ORM + PostgreSQL)
- [HIPPO/dev](https://github.com/xchem/HIPPO/tree/dev): Development branch

</details>


<details>

<summary> Information for HIPPO developers </summary>

### Developer installation

To develop on HIPPO please fork this repository and then install locally:

```bash
git clone https://github.com/YOUR_USER/HIPPO
cd HIPPO
uv sync --frozen
```

Or with pip:

```bash
pip install -e .
```

### Releases

HIPPO is automatically released to [PyPI](https://pypi.org/project/xchem-hippo/) as
`xchem-hippo` via a Github Action off the using the
[release](https://github.com/xchem/HIPPO/actions/workflows/release.yaml) workflow.

When you want to make an official release go to the [Releases](https://github.com/xchem/HIPPO/releases) page
and then click the **Draft a new release** button. Remember to familiarise yourself
with the xchem release process on the trunk-based-development Wiki
[Creating releases](https://github.com/xchem/trunk-based-development/wiki/Creating-releases)
page.

### Code style

HIPPO is linted using [black](https://pypi.org/project/black/) and commits are
automatically linted using the
[lint](https://github.com/xchem/HIPPO/actions/workflows/lint.yaml) workflow.
The use of [pre-commit](https://pre-commit.com/) is encouraged for local development
to automatically run the linting at git commit time:

```bash
pip install pre-commit
pre-commit install
```

### Documentation

Documentation is automatically built off the
[HIPPO/main](https://github.com/xchem/HIPPO/tree/main) branch using readthedocs.
For local building using sphinx:

```bash
cd docs
make html
```

To check API reference coverage use [docstr-coverage](https://pypi.org/project/docstr-coverage/)

```bash
pip install docstr-coverage
docstr-coverage hippo
```

### Tests

Tests require a running PostgreSQL database (see Docker setup below). Run with:

```bash
uv run pytest
```

N.B. the numbered tests, e.g. `test_01_fragalysis_download.py` need to run in sequential order to set up the database. Configure `tests/config.py` to point at your database. The tests will fail if https://fragalysis.diamond.ac.uk can not provide the protein target's data.

</details>

<details>

<summary> Local development with Docker </summary>

### Setting up the environment

HIPPO uses Docker Compose to run PostgreSQL with the RDKit cartridge locally.

1. Create a `.env` file in the project root with your database connection parameters:

```
DB_NAME=designdb
DB_USER=postgres
DB_PASSWORD=your_password
DB_HOST=database
POSTGRES_PORT=5432
```

2. Build the database container (includes RDKit cartridge compilation — this may take some time):

```bash
cd images/xchem-designdb
docker build -t xchem_designdb:latest .
cd ../..
```

3. Build the application container:

```bash
docker build --no-cache . -t hippo_backend:latest
```

4. Launch services:

```bash
docker compose up
```

To run only the database (connecting from your host):

```bash
docker compose up database
```

5. Access the Jupyter environment at the URL printed in the terminal output.

6. Cleanup:

```bash
docker compose down      # stop services
docker compose down -v   # also wipe database volume
```

### Connecting to a remote deployment

Check port availability:

```bash
nc -zv IP_ADDRESS 5432
```

Success will look something like this:

```
Ncat: Version 7.92 ( https://nmap.org/ncat )
Ncat: Connected to IP_ADDRESS:5432.
Ncat: 0 bytes sent, 0 bytes received in 0.01 seconds.
```

To ssh tunnel to a host which has the correct exposed port and forward the correct port:

```bash
ssh -L 5432:IP_ADDRESS:5432 USER@GATEWAY_HOST
```

To test your connection (from your local machine)

```
pg_isready -h localhost -p 5432
```

To list available databases with `psql`

```bash
psql -h localhost -U USER -p 5432 -l
```

To connect to a specific database with `psql`

```bash
psql -h localhost -U USER -p 5432 -n DATABASE
```

</details>
