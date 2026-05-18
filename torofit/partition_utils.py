import numpy as np
import matplotlib.pyplot as plt
from .Fourier_coeffs import *

# def find_partition_for_cocycle_via_statistical_scores(tc, cocycle_idx, thresh=None, std_thresh=0.01):
#     harm_reps = tc._harm_reps
#     from dreimac.utils import CohomologyUtils
#     thresh= tc._rips_threshold if thresh is None else thresh
#     theta = CohomologyUtils.one_cocycle_to_matrix(dist_mat=tc._dist_land_land, 
#                                                 threshold=thresh,
#                                                 lookup_table=tc._cns_lookup_table,
#                                                 n_points = tc._n_landmarks,
#                                                 cocycle_as_vector=harm_reps[cocycle_idx])
#     bad_idx = []
#     good_idx = []
#     #via standard deviation
#     for i in range(len(theta)):
#         vector = theta[i,:]
#         nonzero_vector = vector[vector!=0]
#         if len(nonzero_vector) == 0:
#             bad_idx.append(i)
#         else:
#             diff = np.std(nonzero_vector)
#             if diff <= std_thresh:
#                 bad_idx.append(i)
#             else:
#                 good_idx.append(i)
#     return good_idx, bad_idx


def total_variation_scores(idx, tc, tc_theta, K=None):
    if np.asarray(idx).ndim == 0:
        idx = np.asarray([idx])
    dim = len(idx)
    theta = tc_theta[idx]
    landmarks = tc._idx_land
    dist_land_data = tc._dist_land_data
    r_cover = tc._r_cover
    tvs = []
    for i, landmark in enumerate(landmarks):

        if K is None:
            nbhd = np.where(dist_land_data[i,:] < 2*r_cover )[0] #This is the neighborhood of this landmark
        else:
            nbhd = np.argsort(dist_land_data[i, :])[:K] # K nearest neighbors of this landmark
        tv = 1/dim * np.linalg.norm(np.arccos(np.cos(theta[:, nbhd] - theta[:, landmark][:, None])), axis=0).mean() #This is the total variation score for this landmark
        tvs.append(tv)
    tvs = np.asarray(tvs)
    return tvs

def find_partition_for_cocycle_via_tvs_scores(tc, 
                                              tc_theta, 
                                              cocycle_idx, 
                                              tvs_thresh=0.5):
    scores = total_variation_scores(cocycle_idx, tc, tc_theta)
    thresh = tvs_thresh

    good_idx, bad_idx = np.where(scores > thresh)[0], np.where(scores <= thresh)[0]

    return good_idx, bad_idx

def find_fits(data,
            tc,
            tc_theta,
            cocycle_partitions=None,
            Ks=None, 
            lams=None,
            partition_thresh=None):
    
    if cocycle_partitions is None:
        cocycle_partitions = [[i] for i in range(tc._betti[1])]
    if partition_thresh is None:
        partition_thresh = 0.5*np.ones(len(cocycle_partitions))
    if Ks is None:
        Ks = 5*np.ones(len(cocycle_partitions))
    if lams is None:
        lams = np.ones(len(cocycle_partitions))
    

    fits = []
    for part_id, cocycle_ids in enumerate(cocycle_partitions):
        print(f'partition{cocycle_ids}')
        good_points = np.asarray([]).astype(np.int64)
        bad_points = np.asarray([]).astype(np.int64)
        for i in cocycle_ids:
            
            good_idx, bad_idx = find_partition_for_cocycle_via_tvs_scores(
            tc, 
            tc_theta,
            cocycle_idx=i, 
            tvs_thresh=partition_thresh[i]
        )
        
            landmarks = tc._idx_land
            good_landmarks = landmarks[good_idx]
            mems_for_good = tc._dist_land_data[good_idx, :]
            mems_for_bad = tc._dist_land_data[bad_idx, :]
            
            good_points_cocycle = np.where(np.any(mems_for_good < tc._r_cover, axis=0))[0]
            bad_points_cocycle = np.where(np.any(mems_for_bad < tc._r_cover, axis=0))[0]
            
            good_points = np.append(good_points, good_points_cocycle) #Need union here
            bad_points = np.append(bad_points, bad_points_cocycle)

        coords = np.column_stack([tc_theta[i][good_points] for i in cocycle_ids])
        V = fit(data[good_points], coords, K = Ks[part_id], lam=lams[part_id])
        fits.append(V)
    return fits