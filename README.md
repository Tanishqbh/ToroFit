# ToroFit

[![PyPI version](https://img.shields.io/pypi/v/torofit.svg)](https://pypi.org/project/torofit/)
[![Python versions](https://img.shields.io/pypi/pyversions/torofit.svg)](https://pypi.org/project/torofit/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

ToroFit fits smooth Fourier series on toroidal manifolds $\mathbb{T}^d$ to point-cloud data, using persistent-cohomology circular coordinates from [DREiMaC](https://github.com/Tanishqbh/DREiMac). It is useful for parametrising cyclic structure in any dataset that lives on a (possibly noisy, possibly higher-dimensional) torus — closed curves in $\mathbb{R}^n$, periodic trajectories, branched cycles, or product structures like $S^1 \times S^1$.

> ⚠️ **Alpha software (0.1.0a0).** The API surface — function names, argument order, output shapes — may change in any 0.x release. Pin exact versions if you depend on it.

## Installation

```bash
pip install torofit
```

This pulls a custom branch of DREiMaC (the `partition` branch) automatically.

<details>
<summary>Notes on installation</summary>

ToroFit currently depends on a fork of DREiMaC ([Tanishqbh/DREiMac@partition](https://github.com/Tanishqbh/DREiMac/tree/partition)) for cocycle-partition support that hasn't been upstreamed yet. The fork is declared as a direct URL in `pyproject.toml`, so `pip` resolves it automatically.

Some environments (corporate proxies, conda-forge, restricted CI runners) refuse to install packages from direct git URLs. If `pip install torofit` fails with an error about disallowed direct references, install DREiMaC explicitly first:

```bash
pip install git+https://github.com/Tanishqbh/DREiMac@partition
pip install torofit
```

Once the partition support lands in mainline DREiMaC, this step will disappear.
</details>

## Quickstart

Fit a Fourier series to a noisy closed curve in $\mathbb{R}^3$:

```python
import numpy as np
from torofit import Fourier_coeffs

# Noisy curve in R^3 parameterised by t in [0, 2π).
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

For toroidal coordinates of dimension $d > 1$, pass an `(M, d)` array of angles instead of a 1D `t`. A length-`M` `weights` argument to `fit(...)` enables per-sample weighting in the underlying weighted ridge solve.

End-to-end examples (point cloud → DREiMaC circular coordinates → ToroFit) live in [`docs/notebooks/`](docs/notebooks/).

## Documentation

Full documentation is at [torofit.readthedocs.io](https://torofit.readthedocs.io).

To build locally:

```bash
pip install 'torofit[dev]'
cd docs && make html
```

## Development

```bash
git clone https://github.com/Tanishqbh/ToroFit
cd ToroFit
pip install -e '.[dev]'
pytest
```

The pinned development environment used by the maintainers lives in [`requirements-lock.txt`](requirements-lock.txt). It is not the install path for users — run `pip install torofit` instead.

## Citing

If you use ToroFit in your research, please cite the repository for now; a paper is in preparation. The [`CITATION.cff`](CITATION.cff) file in this repo lets GitHub render a "Cite this repository" widget automatically.

## AI assistance

See [`AI_DISCLOSURE.md`](AI_DISCLOSURE.md) for a statement on the use of AI tools during development of this package.

## License

[MIT](LICENSE) © 2026 Tanishq Bhatia, Jose Perea.
