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
    


