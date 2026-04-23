# ToroFit

A Python package for Topological Data Analysis (TDA).

## Installation

```bash
pip install torofit
```

Or install from source:

```bash
git clone https://github.com/Tanishqbh/ToroFit.git
cd ToroFit
pip install -e .
```

## Quick Start

```python
import torofit

info = torofit.version_info()
print(info)
# {'name': 'torofit', 'version': '0.1.0'}
```

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
