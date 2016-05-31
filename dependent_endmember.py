import burnman
import numpy as np

def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    """

    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out

def indicize(inp_arr):
    arr=np.zeros([len(inp_arr), 1+max(inp_arr.flat)])
    for i, mbr in enumerate(inp_arr):
        arr[i][mbr[0]] = 1
        arr[i][mbr[1]] = 1
    return arr

def dependent_endmembers(formulae):
    solution_formulae, n_sites, sites, n_occupancies, endmember_occupancies, site_multiplicities = burnman.processchemistry.process_solution_chemistry(formulae)

    # First, we redefine our endmembers using unique site occupancy indices
    # We can then find the set of unique occupancies for each site
    i=0
    n_sites = 0
    indices = []
    site_indices = []
    
    for site in sites:
        site_indices.append([])
        site_occupancies = map(tuple, endmember_occupancies[:, i:i+len(site)])
        set_site_occupancies = map(np.array, set(site_occupancies))
        
        for site_occupancy in site_occupancies:
            site_indices[-1].append([i+n_sites for i, set_occupancy in enumerate(set_site_occupancies) if np.array_equal(set_occupancy, site_occupancy)][0])
        indices.append(np.arange(n_sites,n_sites+len(set_site_occupancies)))
        
        i += len(site)
        n_sites += len(set_site_occupancies)

    given_members = indicize(np.array(zip(*site_indices)))

    
    # All the potential endmembers can be created from permutations of the unique site occupancies
    # Note that the solution may not allow all of these endmembers; for example, if there
    # is an endmember with unique occupancies on two sites,
    # (e.g. Ca2+ and Fe3+ in the X, Y sites in garnet)
    # then it will not be possible to exchange Ca2+ for Mg2+ on the X site and retain Fe3+ on the Y site.
    all_members=indicize(cartesian(indices))

    # Now we return only the endmembers which are not contained within the user-provided independent set
    a=np.concatenate((given_members, all_members), axis=0)
    b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    _, idx, cnts = np.unique(b, return_index=True, return_counts=True)
    potential_dependents = a[[i for i, cnt in zip(*[idx, cnts]) if cnt==1]]
    
    
    # We only want the dependent endmembers which can be described by a
    # linear combination of the independent set
    tol = 1e-12
    dependent_endmembers = []
    for mbr in potential_dependents:
        a,resid,rank,s = np.linalg.lstsq(given_members.T, mbr)
        if resid[0] < tol:
            a.real[abs(a.real) < tol] = 0.0
            dependent_endmembers.append(a)

    return np.array(dependent_endmembers)


formulae = ['[Mg]3[Al]2Si3O12',
            '[Mg]3[Mg1/2Si1/2]2Si3O12',
            '[Mg]3[Fe]2Si3O12',
            '[Fe]3[Al]2Si3O12',
            '[Ca]3[Fe1/2Si1/2]2Si3O12']

print dependent_endmembers(formulae)


this = [[0., 1.], [0., 1.]]

for i, (j, k) in enumerate(this):
    print i, j, k
