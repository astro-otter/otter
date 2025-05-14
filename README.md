# OTTER API
### **O**pen mul**T**iwavelength **T**ransient **E**vent **R**epository

A Python API for the OTTER.

[actions-badge]:            https://github.com/astro-otter/otter/workflows/CI/badge.svg
[actions-link]:             https://github.com/astro-otter/otter/actions
[black-badge]:              https://img.shields.io/badge/code%20style-black-000000.svg
[black-link]:               https://github.com/psf/black
[conda-badge]:              https://img.shields.io/conda/vn/conda-forge/hepfile
[conda-link]:               https://github.com/conda-forge/hepfile-feedstock
[github-discussions-badge]: https://img.shields.io/static/v1?label=Discussions&message=Ask&color=blue&logo=github
[github-discussions-link]:  https://github.com/mattbellis/hepfile/discussions
[gitter-badge]:             https://badges.gitter.im/https://github.com/mattbellis/hepfile/community.svg
[gitter-link]:              https://gitter.im/https://github.com/mattbellis/hepfile/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge
[pypi-link]:                https://pypi.org/project/astro-otter/
[pypi-platforms]:           https://img.shields.io/pypi/pyversions/astro-otter
[pypi-version]:             https://badge.fury.io/py/astro-otter.svg
[rtd-badge]:                https://readthedocs.org/projects/otter/badge/?version=latest
[rtd-link]:                 https://otter.readthedocs.io/en/latest/?badge=latest
[sk-badge]:                 https://scikit-hep.org/assets/images/Scikit--HEP-Project-blue.svg
[ruff-badge]:               https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/v2.json
[ruff-link]:                https://github.com/astral-sh/ruff
[codecov-badge]:            https://codecov.io/gh/astro-otter/otter/graph/badge.svg?token=BtCerOdTc0
[codecov-link]:             https://codecov.io/gh/astro-otter/otter

[![Documentation Status](https://readthedocs.org/projects/astro-otter/badge/?version=latest)](https://astro-otter.readthedocs.io/en/latest/?badge=latest)
[![Actions Status][actions-badge]][actions-link]
[![PyPI version][pypi-version]][pypi-link]
[![Linting: Ruff][ruff-badge]][ruff-link]
[![codecov][codecov-badge]][codecov-link]

## Installation
To install the OTTER API use
```
python3 -m pip install astro-otter
```

## Developer Instructions
1. Set the `OTTER_ROOT` environment variable
   ```
   export OTTER_ROOT=/path/to/where/to/clone
   ```
2. Clone the relevant repos:
   ```
   git clone https://github.com/astro-otter/otter.git $OTTER_ROOT/otter
   git clone https://github.com/astro-otter/otterdb.git $OTTER_ROOT/otterdb
   git clone https://github.com/astro-otter/otter-web.git $OTTER_ROOT/otter-web
   ```
3. Install the NASA ADS Python API by following the instructions at https://ads.readthedocs.io/en/latest/#getting-started
4. Install otter, the API for this database. From
   the root directory where you installed these repos:
   ```
   cd $OTTER_ROOT/otter
   python -m pip install -e .
   ```
5. Process the data to build the local copy of the database. Follow the instructions in the otterdb repo README.
6. Easily access the data using the Otter code! In python:
  ```
  import os
  from otter import Otter
  otter = Otter(os.path.join(os.environ['OTTER_ROOT'], 'otterdb', '.otter'))
  res = otter.query(names='AT2018hyz')
  print(res)
  ```

## Installation from Source
To install the OTTER API from the source code use
```
git clone https://github.com/astro-otter/otter.git
cd otter
python -m pip install .
```
This will be changed into the more convenient `python -m pip install astro-otter` at a later date!

For developers, please also enable the pre-commit hooks using
```
pre-commit install
```

## Repo Organization
| Directory | Contents |
|------------|------------|
| `src/otter` | A pip installable API for interfacing with the OTTER database|
| `scripts` | The pipeline scripts for converting unprocessed data into the OTTER JSON format|
| `docs` | Documentation for the OTTER API |
| `test` | Some Unit tests for the source code |
