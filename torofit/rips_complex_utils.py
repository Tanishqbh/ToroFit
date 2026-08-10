import gudhi
from scipy.spatial.distance import cdist
import numpy as np
from scipy.sparse import csr_array


def get_simplex_tree(points, radius, max_dim=2):
    '''
    Args:
        points: coordinates of point cloud
        radius: radius for rips complex
    Returns:
        simplex_tree: simplex_tree of the rips complex at scale radius
    
    '''
    distance_matrix= cdist(points, points)

    rips_complex = gudhi.RipsComplex(distance_matrix=distance_matrix, max_edge_length=radius)

    simplex_tree = rips_complex.create_simplex_tree(max_dimension=max_dim)
    return simplex_tree

def get_d1_for_rips(points, radius):
    '''
    Args:
        points: coordinates of point cloud
        radius: radius for rips complex
    Returns:
        d_1: coboundary matrix in dimension 1
    '''

    simplex_tree = get_simplex_tree(points, radius)
    simplices_2 = [s for s,_ in simplex_tree.get_simplices() if len(s) == 3]
    edges = [s for s,_ in simplex_tree.get_simplices() if len(s) == 2]
    edge_map = {tuple(edge) : i for i, edge in enumerate(edges)}
    data = []
    row_ind = []
    col_ind = []
    for i, simplex in enumerate(simplices_2):
        for j, (tau, _) in enumerate(simplex_tree.get_boundaries(simplex)):
            data.append((-1)**j)
            row_ind.append(i)
            col_ind.append(edge_map[tuple(tau)])
    sparse_d_1 = csr_array((data, (row_ind, col_ind)), shape=(len(simplices_2), len(edges)))
    d_1 = sparse_d_1.todense()
    return d_1
    
        

def get_cup_prod_2_dim(theta1, theta2, points, radius):
    '''
    Args:
        theta1: cocycle 1 
        theta2: cocycle 2
        points: coordinates of point cloud
        radius: radius for rips complex
    Returns:
        cupvec: vector representing theta1 \smile theta 2
    '''
    tree = get_simplex_tree(points, radius)
    simplices_2 = [s for s,_ in tree.get_simplices() if len(s) == 3]
    cupvec = np.zeros(len(simplices_2))
    for i, simplex in enumerate(simplices_2):
        boundary = list(tree.get_boundaries(simplex))
        edge1, edge2 = boundary[0][0], boundary[2][0]
        cupvec[i] = theta1[*edge1] * theta2[*edge2]
    return cupvec

def get_nonzero_cup_products_dim_1(tc, harmonic_class_idx, verbose=False):
    '''
    Args:
        tc: ToroFit object containing the data and parameters
        harmonic_class_idx: list of indices of harmonic classes to consider
        points: coordinates of point cloud
        radius: radius for rips complex
    Returns:
        nonzero_cup_products: list of tuples representing nonzero cup products between 1-dim cohomology classes
    '''
    from dreimac.utils import CohomologyUtils
    harm_reps = tc._harm_reps
    cocycle_mats = []
    for cocycle_idx in harmonic_class_idx:
        theta = CohomologyUtils.one_cocycle_to_matrix(dist_mat=tc._dist_land_land, 
                                                    threshold=tc._rips_threshold,
                                                    lookup_table=tc._cns_lookup_table,
                                                    n_points = tc._n_landmarks,
                                                    cocycle_as_vector=harm_reps[cocycle_idx])
        cocycle_mats.append(theta)
    from itertools import combinations
    cocycle_nonzero_cup_products= []
    points = tc._X[tc._idx_land]; radius= tc._rips_threshold

    d1 = get_d1_for_rips(points, radius=radius)
    
    if verbose:
        print(f"Computing rank for coboundary matrix d1 matrix of size {d1.shape} for Rips complex...")
    rank_d1 = np.linalg.matrix_rank(d1)
    for a,b in combinations(range(len(cocycle_mats)), 2):
        cupvec = get_cup_prod_2_dim(cocycle_mats[a], cocycle_mats[b], points=points, radius=radius)
        rank_cupped = np.linalg.matrix_rank(np.hstack([d1, cupvec.reshape(-1,1)]))
        if rank_d1 < rank_cupped:
            cocycle_nonzero_cup_products.append((a,b))
    return cocycle_nonzero_cup_products 


