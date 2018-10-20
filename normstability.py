import louvainLNL
import louvainLCL
import numpy as np
import numpy.linalg as LA
import scipy 
import scipy.sparse
import scipy.spatial.distance as ssd
#import matplotlib.pyplot as plt
import csv
import pandas as pd

""" Keep this version just in case - more general.  normstabilitymsm.py is
just for MSMBuilder"""

time = 1 #this will be an input parameter 
precision = 1e-9 #also input

teleport_tau = 0.15
nnodes = adjacency.shape[0]

#may not need this if guaranteed system is ergodic
if directed == True:
    dout = np.sum(adjacency, axis=1)
    dangling = (dout==0)
    dout(dangling) = 1
    Dout = np.diag(dout)

    M = (1 - teleport_tau)*np.dot(np.linalg.inv(Dout), adjacency) #may be better to use linear solve here
    #teleportation according to arXiv:0812.1770
    M = M + np.dot(np.diag(teleport_tau + dangling*(1 - teleport_tau)), 
                       np.ones((nnodes,nnodes))/nnodes)

#if not then following code is fine for directed
dout = np.sum(adjacency, axis=1)
Dout = np.diag(dout)

M = np.dot(np.linalg.inv(Dout), adjacency) #may be better to use linear solve here
u,v = np.linalg.eig(M.T) #must be M_transpose if originally defined as M_ij : i --> j

lambda_ = np.max(u)
pi_norm = v[:, u == np.max(u)] #extract column corresponding to eigenvalue = 1
pi_norm = np.abs(statdist) #make sure all values are positive

pi = pi_norm/np.sum(pi_norm)

#compute exponential transition matrix
solution = np.diag(pi).dot(scipy.linalg.expm(time*(M - np.eye(len(M), len(M)))))

#if MSMBuilder generates the transition rate matrix (K) then start here:
solution = np.diag(pi).dot(scipy.linalg.expm(time*K))

#symmetrize
solution = (solution+solution.T)/2


#prune weights that are too small according to precision parameter
solution = (np.max(solution)*precision)*np.round(solution/((np.max(solution)*precision)))
graph = find(solution)


(stability, number_of_comms, community_id, VI) = _full_stability(graph, louvain_runs, precision, calcVI)




def _full_stability(graph, louvain_runs, precision, calcVI):

    louvain_ensemble = []
    stability_partition_list = []
    number_of_comms_partition_list = []

    for i in range(louvain_runs):
        (stability, number_of_comms, community_id) = louvainLNL.stability(graph, 1, precision)
        
        louvain_ensemble.append(community_id)
        stability_partition_list.append(stability)
        number_of_comms_partition_list.append(number_of_comms)
        

    stability = max(stability_partition_list)
    index = stability_partition_list.index(max(stability_partition_list))
    number_of_comms = number_of_comms_partition_list[index]
    community_id = louvain_ensemble[index]

    if calcVI:
        VI = varinfo(louvain_ensemble)  # Added here
        return (stability, number_of_comms, community_id, VI) #will also need to return louvain ensemble
    
    else:
        VI = [] #seems a bit hacky, find a way to return different numbers of variables 
        return (stability, number_of_comms, community_id, VI)



def varinfo(louvain_ensemble):

    louvain_ensemble = np.array(louvain_ensemble) #convert to numpy array

    number_of_partitions = louvain_ensemble.shape[0]
    n = louvain_ensemble.shape[1]
    VI_mat = np.zeros(number_of_partitions)
    VI = 0


    # If all the partitions are identifcal, VI = 0 and there is no need to do the
    # rest of the calculations which are computationally expensive

    if np.all( louvain_ensemble == np.tile(louvain_ensemble[0,:], (number_of_partitions, 1))):
        return 0 # np.zeros((number_of_partitions, number_of_partitions))) - this refers to VI matrix

    # select only partitions that are different
    # need to copy this into new array to stop numpy converting tuples to arrays
    # otherwise 'np.unique' flattens result into single array rather than returning
    # unique rows
    louvain_ensemble = [tuple(row) for row in louvain_ensemble]
    temp_louvain_ensemble = np.empty((number_of_partitions,), dtype=object)
    # enumerate doesn't work here, creates arrays
    for i in range(len(louvain_ensemble)):
        temp_louvain_ensemble[i] = louvain_ensemble[i]

    (louvain_ensemble, indices)  = np.unique(temp_louvain_ensemble, return_inverse = True)
    number_of_partitions = len(louvain_ensemble)
    VI_mat = np.zeros((number_of_partitions, number_of_partitions))

    VI_tot = 0
    nodes = np.linspace(0,n-1,n)
    ones = np.ones(n)

    # need to look into parallelizing loops in numpy
    for i in range(number_of_partitions):
        A_1 = scipy.sparse.csr_matrix((ones , (louvain_ensemble[i], nodes)))
        n_1_all = np.sum(A_1.toarray(), axis=1)
        
        for j in range(i):
            A_2 = scipy.sparse.csr_matrix((ones,(nodes, louvain_ensemble[j])))
            n_2_all = np.sum(A_2.toarray(), axis=0)#check
            n_12_all = np.dot(A_1, A_2)

            (rows, cols, n_12) = find(n_12_all.toarray())

            n_1 = n_1_all[np.array(rows).astype(int)]
            n_2 = n_2_all[np.array(cols).astype(int)]

            VI = np.sum(n_12*np.log((n_12**2)/(n_1 * n_2))) # all element-wise operations
            VI = -1/(n*np.log(n))*VI
            VI_mat[i,j] = VI
            VI_tot = VI_tot+VI
        
    
    VI_mat_full = np.zeros((number_of_partitions, len(indices)))

    for i in range(number_of_partitions):
        VI_mat_full[i] = VI_mat[i, np.array(indices)]
    
    VI_mat_full = VI_mat_full[np.array(indices)]

    VI_mat = VI_mat_full + np.transpose(VI_mat_full)

    VI = np.mean(ssd.squareform(VI_mat))

    return VI # decide on VI_mat









# helper function - replicates Matlab's find function.  Takes ndarray as argument.
def find(matrix):

    unweighted_partition_list = list(zip(np.where(matrix > 0))) #gives [row, col] without values

    edge_partition_list = np.where(matrix > 0)
    partition_list_of_indices = list(zip(edge_partition_list[0], edge_partition_list[1])) #i.e. [(0, 1), (1, 0), (1, 2), (1, 3), (2, 1), (2, 3), (3, 1), (3, 2)]

    values = [matrix[element[0]][element[1]] for element in partition_list_of_indices] #extracts weighting values using indices above

    graph = np.array([unweighted_partition_list[0][0], unweighted_partition_list[1][0], values], dtype=np.float64) #bit ugly, should fix

    return graph

















