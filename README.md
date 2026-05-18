# ToroFit

## Description

This package implements the toroidal fitting algorithm using dreimac.

## Installation

This package depends on a custom branch of dreimac, namely the [partition branch](https://github.com/Tanishqbh/DREiMac/tree/partition).
```
pip install git+https://github.com/Tanishqbh/DREiMac@partition
pip install -e .
```

## Quick Start

Fit a Fourier series to a noisy closed curve in $\mathbb{R}^3$:

```python
import numpy as np
from torofit import Fourier_coeffs

# Sample a noisy curve in R^3 parameterised by t in [0, 2π).
t = np.linspace(0, 2 * np.pi, 200, endpoint=False)
X = np.column_stack([
    np.cos(t),
    np.sin(t),
    0.5 * np.cos(2 * t),
]) + 0.05 * np.random.randn(200, 3)

# Fit a Fourier series with bandwidth K and Dirichlet-energy regulariser lam.
K, lam = 5, 1.0
V = Fourier_coeffs.fit(X, t, K=K, lam=lam)

# Evaluate the fitted curve at new angles.
t_new = np.linspace(0, 2 * np.pi, 1000).reshape(-1, 1)
fitted = Fourier_coeffs.Toroidal_Evaluate(V, t_new, K=K)
```

For toroidal coordinates `Y` of dimension `d > 1`, pass an `(M, d)` array instead of a 1D `t`. A length-`M` `weights` argument enables per-sample weighting in the underlying weighted least-squares solve.

## Documentation

Full documentation is available at [torofit.readthedocs.io](https://torofit.readthedocs.io).

To build the docs locally:

```bash
pip install sphinx
cd docs
make html
```

## Development

Run the tests:

```bash
pip install pytest
pytest
```

## License

MIT
