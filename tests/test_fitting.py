from torofit import Fourier_coeffs
import numpy as np
from dreimac import ToroidalCoords
def sample_noisy_torus(R, r, num_points, noise_std):

    u = np.random.uniform(0, 2*np.pi, num_points)

    v = np.random.uniform(0, 2*np.pi, num_points)

    x = (R + r * np.cos(u)) * np.cos(v) + np.random.normal(scale=noise_std, size=num_points)

    y = (R + r * np.cos(u)) * np.sin(v) + np.random.normal(scale=noise_std, size=num_points)

    z = r * np.sin(u) + np.random.normal(scale=noise_std, size=num_points)

    return np.array([x, y, z]).T

class TestFitting:
    def test_fit_computation(self):
        data = sample_noisy_torus(5.0,1.0,10000,0.0)
        tc = ToroidalCoords(data,n_landmarks=2000)
        tc_theta = tc.get_coordinates(perc=0.0,cocycle_idxs=[0,1]).T
        
        X = complex(1,0)*data
        Y = tc_theta
        K = 4
        lam = 1.0
        V=Fourier_coeffs.fit(X,Y,K,lam)
        
class TestComplexLinAlg:
    def test_least_squares(self):
        i = complex(0,1)
        A = np.asarray([[3+i, 0, 6, 2],
                        [2, 1, 0, 0],
                        [1, 0, 1+2*i, 0],
                        [1, 1, 1, 1]])
        b = np.asarray([[1],
                        [i],
                        [1],
                        [1]])
        x = Fourier_coeffs.solve_least_squares_(A,b)
        assert np.allclose(A@x,b)


class TestWeightedFit:
    def _toy_problem(self, seed=0):
        from torofit.Fourier_coeffs import Exponentiated_coordinates, Energy_weights
        rng = np.random.default_rng(seed)
        M, d, K, N = 200, 2, 3, 4
        Y = rng.uniform(0, 2 * np.pi, size=(M, d))
        X = rng.normal(size=(M, N))
        E = Exponentiated_coordinates(Y, K, d)
        Q = np.diag(Energy_weights(K, d) ** 2)
        return X, Y, K, E, Q, M

    def test_none_matches_uniform_ones(self):
        from torofit import Fourier_coeffs
        X, Y, K, _, _, M = self._toy_problem()
        lam = 0.5
        V_none = Fourier_coeffs.fit(X, Y, K, lam)
        V_ones = Fourier_coeffs.fit(X, Y, K, lam, weights=np.ones(M))
        assert np.allclose(V_none, V_ones)

    def test_normal_equations(self):
        from torofit import Fourier_coeffs
        X, Y, K, E, Q, M = self._toy_problem(seed=1)
        lam = 0.5
        w = np.random.default_rng(2).uniform(0.1, 2.0, size=M)
        V = Fourier_coeffs.fit(X, Y, K, lam, weights=w)
        W = np.diag(w)
        lhs = (E.conj().T @ W @ E + lam * Q) @ V
        rhs = E.conj().T @ W @ X
        assert np.allclose(lhs, rhs)

    def test_scale_invariance(self):
        # Scaling all weights by c is equivalent to scaling lam by 1/c.
        from torofit import Fourier_coeffs
        X, Y, K, _, _, M = self._toy_problem(seed=3)
        lam, c = 0.5, 2.0
        V_scaled_w = Fourier_coeffs.fit(X, Y, K, lam, weights=c * np.ones(M))
        V_scaled_lam = Fourier_coeffs.fit(X, Y, K, lam / c)
        assert np.allclose(V_scaled_w, V_scaled_lam)