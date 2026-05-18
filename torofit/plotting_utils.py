import numpy as np
from itertools import product
from .Fourier_coeffs import *
def plot_fitted_torus_via_coeff_matrix(V, K, samples=100):
    '''
    Parameters
    ----------
    V : np.ndarray (2K+1)**d X N, coefficient matrix
    K : int, bandwidth 

    Return
    ------
    image : vertices of mesh
    triangles : triangles of the mesh (Is,Js,Ks) to be used for plotting via plotly
    '''
    Is, Js, Ks = [], [], []
    N = samples
    theta = np.linspace(0, 2*np.pi, num=N)
    domain = np.asarray(list(product(theta, repeat=2))).T
    image = Toroidal_Evaluate(V,domain.T,K).real
    for i in np.arange(N):
        for j in np.arange(N):
            coord = i + N*j
            Is = Is + [coord,coord]
            Js = Js + [(coord + 1) % (N**2), (coord + N) % (N**2)]
            Ks = Ks + [(coord + 1 + N) % (N**2), (coord + 1 + N) % (N**2)]
    triangles = np.column_stack([Is, Js, Ks])

    return image, triangles