Installation
============

Requirements
------------

* Python 3.9, 3.10, or 3.11 (Python 3.12+ is not supported)
* pip

ToroFit depends on a custom branch of `DREiMac
<https://github.com/Tanishqbh/DREiMac/tree/partition>`_ (the
``partition`` branch), which provides the toroidal coordinates and
support-partition utilities used by the fitting pipeline. DREiMac
itself requires Python ``<3.12``, which is what bounds ToroFit's
supported Python range.

Conda environment
-----------------

A ready-made ``environment.yml`` is provided at the repository root.
It pins the supported Python range and installs all runtime
dependencies via conda-forge plus the DREiMac partition branch and
an editable ToroFit install via pip:

.. code-block:: bash

   conda env create -f environment.yml
   conda activate torofit

Install from source
-------------------

.. code-block:: bash

   pip install git+https://github.com/Tanishqbh/DREiMac@partition
   git clone https://github.com/Tanishqbh/ToroFit.git
   cd ToroFit
   pip install -e .
