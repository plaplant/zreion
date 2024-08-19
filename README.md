# zreion

`zreion` is a way of quickly computing a "redshift of reionization" field using
a semi-numeric method developed in [Battaglia et
al. (2013)](https://ui.adsabs.harvard.edu/abs/2013ApJ...776...81B/abstract). The
method assumes the redshift of reionization for a particular point in a
cosmological volume is a biased tracer of the matter field, and can be written
as a parameterized bias function. For a full derivation and comparison with
simulations, see the linked paper.

## Installation

Installing the package can be performed by checking the repo out, and then running:

```
FC="ifort" pip install .
```

Dependencies should be handled automatically by `pip` if not already installed.

## Dependencies

* setuptools >= 38.3
* numpy
* pyfftw

## Tests

The repo includes a test suite, which can be invoked by running `pytest` in the
top-level of the repo. In addition to the package dependencies, it requires
`packaging`, `pytest-cov` and `pytest-cases`. These can be installed using:

```
pip install .[test]
```