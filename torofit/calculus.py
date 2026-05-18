from .Fourier_coeffs import Energy_weights, Toroidal_Evaluate, Exponentiated_coordinates
import numpy as np

def derivative_coeffs(V, K):
    idx = np.arange(len(V)) - K
    return complex(0, 1)*V * idx[:, np.newaxis]

def arclength_func(V, K):
    '''
    Parameters
    ----------

    V : Coeffient matrix for 1d fourier series
    K : Number of modes

    Returns
    -------

    s_from_t : s(t), arclength function
    t_from_s : t(s), inverse of arclength function
    total_length: total length of the curve

    '''
    import numpy as np
    from scipy.integrate import quad
    from scipy.interpolate import interp1d
    
    # Define a simple 2D curve
    
    def r(t):
        res = Toroidal_Evaluate(V, np.asarray([[t]]), K)
        return res.real.reshape((-1,))
    
    # Calculate the derivative
    dV = derivative_coeffs(V, K)
    
    
   
    def r_prime(t):
        et = Exponentiated_coordinates(np.asarray([[t]]), K, 1)
        res = et @ dV
        return res.real.reshape((-1,))

    # Calculate the magnitude of the derivative
    def magnitude_r_prime(t):
        return np.linalg.norm(r_prime(t))

    # Calculate the arc length function (numerical integration)
    def s_function(t, t0=0):
        result, _ = quad(magnitude_r_prime, t0, t)
        return result

    # Generate a lookup table for (t, s)
    t_values = np.linspace(0, 2*np.pi, 100, endpoint=True)
    s_values = [s_function(t) for t in t_values]

    # Create an interpolator for t(s)
    t_from_s = interp1d(s_values, t_values, kind='linear')
    s_from_t = interp1d(t_values, s_values, kind='linear')
    return s_from_t, t_from_s, s_values[-1]