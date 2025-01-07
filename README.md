# README

Routines for making use of superancillary equations and K-D trees to make iterative calculations with thermodynamic models much faster and much more reliable.

Much more information in the docs.

## Docs

Build of docs follow the instructions in the ``doc`` folder

## Run tests

```
mkdir bld
cd bld
cmake .. -DTEQPFLSH_TESTS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release
```
All tests should pass on all platforms

## Development information

### How to do incremental rebuilds with nanobind

See https://nanobind.readthedocs.io/en/latest/packaging.html#step-5-incremental-rebuilds ::

```
pip install nanobind scikit-build-core[pyproject]
pip install --no-build-isolation -ve .
```
