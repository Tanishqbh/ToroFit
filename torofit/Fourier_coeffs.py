import numpy as np
from scipy.linalg import solve_triangular


def solve_least_squares_(A, b):
    """Solve least-squares Ax = b via Cholesky on the normal equations."""
    AH = A.conj().T
    L = np.linalg.cholesky(AH @ A)
    z = solve_triangular(L, AH @ b, lower=True)
    return solve_triangular(L.conj().T, z, lower=False)


def flatten_grid(K, d):
    """All integer frequency multi-indices in [-K, K]^d, shape (F, d)."""
    axes = np.arange(-K, K + 1)
    grids = np.meshgrid(*[axes] * d, indexing='ij')
    return np.column_stack([g.ravel() for g in grids])


def Exponentiated_coordinates(Y, K, d):
    """Fourier feature matrix exp(i L yᵀ) for each row y in Y, shape (M, F)."""
    L = flatten_grid(K, d)
    return np.exp(1j * np.asarray(Y) @ L.T)


def Energy_weights(K, d):
    """L2 norm of each frequency multi-index, shape (F,)."""
    if d == 0:
        return np.ones(1)
    return np.linalg.norm(flatten_grid(K, d), axis=1)


def fit(X, Y, K, lam, weights=None):
    """Fit a T^d Fourier series to data X using toroidal coordinates Y.

    Solves the weighted ridge problem
        argmin_Ω  ||W^{1/2}(Θ_K Ω - X)||_F^2 + λ ||Ω||_Q^2,
    where W = diag(weights). When weights is None, all samples get weight 1.
    """
    if not np.iscomplexobj(X):
        X = X.astype(complex)

    Y = np.asarray(Y)
    if Y.ndim < 2:
        Y = Y.reshape(-1, 1)

    d = Y.shape[1]
    E = Exponentiated_coordinates(Y, K, d)
    r = Energy_weights(K, d)

    if weights is not None:
        w_sqrt = np.sqrt(np.asarray(weights, dtype=float)).reshape(-1, 1)
        E = w_sqrt * E
        X = w_sqrt * X

    F = (2 * K + 1) ** d
    E_prime = np.vstack([E, np.sqrt(lam) * np.diag(r)])
    X_prime = np.vstack([X, np.zeros((F, X.shape[1]))])

    return solve_least_squares_(E_prime, X_prime)


def Toroidal_Evaluate(V, angles, K):
    """Evaluate Fourier series with coefficients V at rows of angles, shape (M, N)."""
    d = angles.shape[1]
    return (Exponentiated_coordinates(angles, K, d) @ V).real
