""" 
For protein stability, we only use combinatorial stability,
with weighted, undirected networks so only parameters required
are:

- Number of Louvain optimisations
- Full or Linearised stability
- Precision
- Parallel -> haven't sorted this out yet


"""

""" Full Stability - uses louvainLNL """

import louvainLNL
import louvainLCL
import numpy as np
import scipy 
import scipy.sparse
import scipy.spatial.distance as ssd
#import matplotlib.pyplot as plt
import csv
import pandas as pd




"""
Need precision and louvain runs here now
"""

def calculate_full_stability(adjacency, time_array, louvain_runs=100, precision=1e-9,  calcVI=True):
    
    
    timesteps = [element[0] for element in enumerate(time_array)]
    
    no_of_nodes = adjacency.shape[0]

    avg_degree = np.sum(adjacency)/no_of_nodes

    # L = D - A
    Laplacian =  np.diag(np.sum(adjacency, axis=0)) - adjacency

    PI = np.diag(np.ones(no_of_nodes)/no_of_nodes)


    stability_array = []
    number_of_comms_array = []
    community_id_array = []
    VI_array = []

    
    #uses louvainLNL
    for i,time in enumerate(time_array): #time is an array of times

        exponential = scipy.linalg.expm(-time*Laplacian/avg_degree)

        solution = np.dot(PI, exponential)

        #need to strip out matrix elements according to provided precision value
        solution = np.max(solution) * precision * np.round(solution/(np.max(solution)*precision))

        graph = find(solution) #find() is a helper function - see below

        (stability, number_of_comms, community_id, VI) = _full_stability(graph, louvain_runs, precision, calcVI)

        stability_array.append(stability)
        number_of_comms_array.append(number_of_comms)
        community_id_array.append(community_id)
        VI_array.append(VI)

        print ('Timestep ' + str(i+1) + ' of ' + str(len(time_array)) + ' completed.')            


    stability_results_frame = pd.DataFrame(
        {
            'Markov time' : time_array,
            'stability' : stability_array,
            'number_of_communities' : number_of_comms_array,
            'community_id' : community_id_array,
            'VI' : VI_array
        },
        index=timesteps,
    )

    return stability_results_frame
    

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





def calculate_linear_stability(adjacency, time_array, louvain_runs=100, precision=1e-9, calcVI=True):
    """ 
    Need adjacency to be a dense array here
    """      


    timesteps = [element[0] for element in enumerate(time_array)]

    graph = find(adjacency) 

    stability_array = []
    number_of_comms_array = []
    community_id_array = []
    VI_array = []

    for i, time in enumerate(time_array): #time is an array of times

        (stability, number_of_comms, community_id, VI) = _linear_stability(graph, time, louvain_runs, precision, calcVI)

        stability_array.append(stability)
        number_of_comms_array.append(number_of_comms)
        community_id_array.append(community_id)
        VI_array.append(VI)

        print ('Timestep ' + str(i+1) + ' of ' + str(len(time_array)) + ' completed.')

    stability_results_frame = pd.DataFrame(
        {
            'Markov time' : time_array,
            'stability' : stability_array,
            'number_of_communities' : number_of_comms_array,
            'community_id' : community_id_array,
            'VI' : VI_array
        },
        index=timesteps,
    )

    return stability_results_frame            


def _linear_stability(graph, time, louvain_runs, precision, calcVI):

    louvain_ensemble = []
    stability_partition_list = []
    number_of_comms_partition_list = []

    for i in range(louvain_runs):
        (stability, number_of_comms, community_id) = louvainLCL.stability(graph, time, precision)
        
        louvain_ensemble.append(community_id)
        stability_partition_list.append(stability)
        number_of_comms_partition_list.append(number_of_comms)
        

    stability = max(stability_partition_list)
    index = stability_partition_list.index(max(stability_partition_list))
    number_of_comms = number_of_comms_partition_list[index]
    community_id = louvain_ensemble[index]

    if calcVI:
        VI = varinfo(louvain_ensemble) # Added here
        return (stability, number_of_comms, community_id, VI)

    else:
        VI = [] #seems a bit hacky, find a way to return different numbers of variables 
        return (stability, number_of_comms, community_id, VI)



"""
#basic plot function
def plot_stability(timesteps, number_of_communities, stability, VI):

    x = np.array(timesteps)
    y1 = np.array(number_of_communities)
    y2 = np.array(stability)
    y3 = np.array(VI)

    ax1 = plt.subplot(211)
    ax1.plot(x, y1, 'r')
    ax1.set_ylabel('Number of communities')

    ax2 = ax1.twinx()
    ax2.plot(x, y2, 'g')
    ax2.set_ylabel('Stability')

    ax3 = plt.subplot(212)
    ax3.plot(x, y3, 'b')
    ax3.set_xlabel('Markov time (s)')
    ax3.set_ylabel('Variation of Information')

    plt.show()

"""    
    

def csv_for_sankey(timesteps=[], filename=''):

    """ Enter an array of the timesteps you wish to use in the Sankey diagram and the filename of the output file """
    
    C = np.array(self.community_id_array)

    partitions = {}
    for n,t in enumerate(timesteps):

        comms = int(max(C[t,:]) + 1) # gives the index to loop over within a timestep (i.e. number of communities - 1)

        d={}
        for x in range(0, comms):
            d["{0}".format(x)] = []


        for i,k in enumerate(C[t,:]):
            for j in d.keys(): 
                if k == int(j):
                    d[j].append(str(i))

        partitions["timestep {0}".format(t)] = d

    #Need now to generate a partition_list of 'values' between the communities in different timesteps
    #Number of timesteps again given by a and b

    partition_list = []
    for n,t in enumerate(timesteps): 
        

        initial_community =  partitions["timestep {0}".format(timesteps[n])] 
        final_community  =  partitions["timestep {0}".format(timesteps[n+1])]

        initial_comms = int(max(C[timesteps[n],:]) + 1)  
        final_comms = int(max(C[timesteps[n+1],:]) + 1)

        for i in range(0, initial_comms):
            for j in range(0, final_comms):
                
                initial_partition_list = initial_community[str(i)]
                final_partition_list = final_community[str(j)]

                value = len(set(initial_partition_list).intersection(final_partition_list))

                each_tuple = ("T {0}: C {1}".format(timesteps[n],i), 
                        "T {0}: C {1}".format(timesteps[n+1], j),
                        value)
                partition_list.append(each_tuple)

        if n+1 == len(timesteps) - 1:
            break

    final_partition_list = [value for value in partition_list if value[2] != 0]

    with open(filename,'w') as out: #change output file here
        csv_out=csv.writer(out)
        csv_out.writerow(['source','target', 'value'])

        for row in final_partition_list:
            csv_out.writerow(row)
     

    

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


