import numpy as np
from .Fourier_coeffs import flatten_grid, Toroidal_Evaluate, Exponentiated_coordinates, Energy_weights


def get_energy(coeff, K, d):
    r = Energy_weights(K, d)
    return np.linalg.norm(r[:, None] * coeff)

def get_fitting_quality(coeff, cc_theta, X, K, d):
    res = Toroidal_Evaluate(coeff, cc_theta, K).real
    quality = np.linalg.norm(res - X)
    return quality


def compute_condition_number(components, wedge_point):
    '''
    Parameters
    -------------------------------------------------
    components: List[np.ndarray], list of component spaces
    wedge_point: np.ndarray, coordinates of the wedge point

    Returns
    -------------------------------------------------
    kappa: float, condition number of the wedge of the component spaces
    
    '''
    from scipy.spatial.distance import cdist
    
    if len(components) < 2:
        raise ValueError("components must contain at least two component spaces")

    wedge_point = np.asarray(wedge_point).reshape(1, -1)
    normalized_components = [np.asarray(comp) for comp in components]

    if any(comp.ndim != 2 for comp in normalized_components):
        raise ValueError("each component space must be a 2D array")

    if any(len(comp) == 0 for comp in normalized_components):
        raise ValueError("component spaces must be non-empty")

    dists_to_wedge_point = [cdist(comp, wedge_point).reshape(-1,) for comp in normalized_components]
  
    blocks = []
    ## We go i,j first
    for i, comp1 in enumerate(normalized_components):
        block = []
        for j, comp2 in enumerate(normalized_components):
            if i == j:
                continue
            else:
                ## Distance between a point in ith component and jth component is simply the distance through the wedge point
                mat = dists_to_wedge_point[i][:, None] + dists_to_wedge_point[j][None, :]
                ambient_dists = cdist(comp1, comp2)
                block.append(np.max(mat - ambient_dists))
        block = np.max(block)
        blocks.append(block)

    return np.max(blocks)

def wedge_defects(components, wedge_point):
    from scipy.spatial.distance import cdist
    import numpy as np

    wedge_point = np.asarray(wedge_point).reshape(1, -1)

    d_to_p = [
        cdist(comp, wedge_point).ravel()
        for comp in components
    ]

    defects = []

    for i in range(len(components)):
        for j in range(i + 1, len(components)):

            ambient = cdist(components[i], components[j])

            through_p = (
                d_to_p[i][:, None]
                + d_to_p[j][None, :]
            )

            defects.append((through_p - ambient).ravel())

    return np.concatenate(defects)

def wedge_defect_profile(
    components,
    wedge_point,
    distance_mode="euclidean",
    geodesic_n_landmarks=400,
    geodesic_k=10,
):
    """
    Parameters
    -------
    components : List[np.ndarray]
        Component point clouds.

    wedge_point : np.ndarray
        Coordinates of the wedge point.

    distance_mode : str, optional
        Distance backend for pairwise distances. Supported values are
        "euclidean" (default, uses scipy cdist) and "geodesic"
        (uses dreimac landmark geodesic distances).

    geodesic_n_landmarks : int, optional
        Number of landmarks used by dreimac geodesic routine.

    geodesic_k : int, optional
        Number of neighbors used by dreimac geodesic routine.

    Returns
    -------
    ambient_scales : np.ndarray
        Ambient cross distances.

    defects : np.ndarray
        Corresponding wedge defects.
    """

    import numpy as np
    from scipy.spatial.distance import cdist

    if distance_mode not in {"euclidean", "geodesic"}:
        raise ValueError("distance_mode must be 'euclidean' or 'geodesic'")

    components = [np.asarray(comp) for comp in components]

    if len(components) < 2:
        raise ValueError("components must contain at least two component spaces")

    if any(comp.ndim != 2 for comp in components):
        raise ValueError("each component space must be a 2D array")

    if any(len(comp) == 0 for comp in components):
        raise ValueError("component spaces must be non-empty")

    wedge_point = np.asarray(wedge_point).reshape(1, -1)

    if any(comp.shape[1] != wedge_point.shape[1] for comp in components):
        raise ValueError("all component spaces and wedge_point must have matching dimensions")

    def _euclidean_dists(comp_a, comp_b):
        return cdist(comp_a, comp_b)

    def _geodesic_dists(comp_a, comp_b):
        try:
            from dreimac import GeometryUtils
        except ImportError as exc:
            raise ImportError(
                "Geodesic mode requires dreimac. Install it to use distance_mode='geodesic'."
            ) from exc

        X = np.vstack([comp_a, comp_b])
        n_a = len(comp_a)
        n_b = len(comp_b)

        dist_mat, _ = GeometryUtils.landmark_geodesic_distance(
            X,
            geodesic_n_landmarks,
            geodesic_k,
        )

        dist_mat = np.asarray(dist_mat)
        if dist_mat.ndim != 2:
            raise ValueError("dreimac geodesic routine returned a non-2D distance matrix")

        n_total = n_a + n_b

        if dist_mat.shape == (n_total, n_total):
            return dist_mat[:n_a, n_a:]

        # dreimac can return landmark-to-point distances with shape
        # (n_landmarks, n_points) or (n_points, n_landmarks).
        if dist_mat.shape[1] == n_total:
            landmark_to_points = dist_mat
        elif dist_mat.shape[0] == n_total:
            landmark_to_points = dist_mat.T
        else:
            raise ValueError(
                "dreimac geodesic distance matrix has unexpected shape for the given inputs"
            )

        A = landmark_to_points[:, :n_a]
        B = landmark_to_points[:, n_a:]

        cross_dists = np.full((n_a, n_b), np.inf, dtype=float)

        for l in range(landmark_to_points.shape[0]):
            cross_dists = np.minimum(cross_dists, A[l][:, None] + B[l][None, :])

        return cross_dists

    dist_fn = _euclidean_dists if distance_mode == "euclidean" else _geodesic_dists

    d_to_p = [
        dist_fn(comp, wedge_point).ravel()
        for comp in components
    ]

    ambient_scales = []
    defects = []

    for i in range(len(components)):
        for j in range(i + 1, len(components)):

            ambient = dist_fn(components[i], components[j])

            through_p = (
                d_to_p[i][:, None]
                + d_to_p[j][None, :]
            )

            defect = through_p - ambient

            ambient_scales.append(ambient.ravel())
            defects.append(defect.ravel())

    return (
        np.concatenate(ambient_scales),
        np.concatenate(defects),
    )
def local_kappa_function(ambient_scales, defects):
    """
    Returns a callable kappa(r).
    """

    import numpy as np

    order = np.argsort(ambient_scales)

    scales = ambient_scales[order]
    defects_sorted = defects[order]

    running_max = np.maximum.accumulate(defects_sorted)

    def kappa(r):

        idx = np.searchsorted(scales, r, side='right') - 1

        if idx < 0:
            return 0.0

        return running_max[idx]

    return kappa