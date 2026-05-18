import numpy as np
from .Fourier_coeffs import flatten_grid, Toroidal_Evaluate, Exponentiated_coordinates, Energy_weights


def get_energy(coeff, K, d):
    r = Energy_weights(K, d)
    return np.linalg.norm(r[:, None] * coeff)

def get_fitting_quality(coeff, cc_theta, X, K, d):
    res = Toroidal_Evaluate(coeff, cc_theta, K).real
    quality = np.linalg.norm(res - X)
    return quality

