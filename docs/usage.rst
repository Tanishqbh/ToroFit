Usage
=====

Fitting a Fourier series to a closed curve
------------------------------------------

The core fitting routine :func:`torofit.Fourier_coeffs.fit` solves the
weighted ridge problem

.. math::

   \Omega^{*} = \arg\min_{\Omega}\;
   \lVert W^{1/2}(\Theta_{K}\Omega - X)\rVert_{F}^{2}
   + \lambda\,\lVert\Omega\rVert_{Q}^{2},

where :math:`\Theta_{K}` is the Fourier feature matrix of bandwidth
:math:`K`, :math:`Q` is the diagonal Dirichlet-energy weighting, and
:math:`W` is an optional diagonal of per-sample weights.

.. code-block:: python

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

Higher-dimensional tori
-----------------------

For :math:`\mathbb{T}^{d}` with :math:`d > 1`, pass an ``(M, d)`` array of
toroidal coordinates as ``Y``. In practice these come from the toroidal
coordinates algorithm of Scoccola et al., available via
:class:`dreimac.ToroidalCoords` on the ``partition`` branch (see the
installation page).

Per-sample weights
------------------

Pass a length-``M`` array via the optional ``weights`` argument to
downweight individual samples — for example by inverse local density,
total-variation score, or distance to a base point:

.. code-block:: python

   V = Fourier_coeffs.fit(X, t, K=K, lam=lam, weights=w)

Omitting ``weights`` (or passing ``None``) reproduces the unweighted fit
exactly.
