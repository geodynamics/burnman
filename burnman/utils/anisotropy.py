import numpy as np

try:  # numpy.block was new in numpy version 1.13.0.
    block = np.block(
        [
            [np.ones((3, 3)), 2.0 * np.ones((3, 3))],
            [2.0 * np.ones((3, 3)), 4.0 * np.ones((3, 3))],
        ]
    )
except:
    block = np.array(
        np.bmat(
            [[[[1.0] * 3] * 3, [[2.0] * 3] * 3], [[[2.0] * 3] * 3, [[4.0] * 3] * 3]]
        )
    )
voigt_compliance_factors = block


def voigt_index_to_ij(m):
    """
    Returns the ij (or kl) indices of the
    stiffness tensor which correspond to those
    of the Voigt notation m (or n).
    """
    if m == 3:
        return 1, 2
    elif m == 4:
        return 0, 2
    elif m == 5:
        return 0, 1
    else:
        return m, m


def voigt_notation_to_stiffness_tensor(voigt_notation):
    """
    Converts a stiffness tensor in Voigt notation (6x6 matrix)
    to the full fourth rank tensor (3x3x3x3 matrix).
    """
    stiffness_tensor = np.zeros([3, 3, 3, 3])
    for m in range(6):
        i, j = voigt_index_to_ij(m)
        for n in range(6):
            k, l = voigt_index_to_ij(n)
            stiffness_tensor[i][j][k][l] = voigt_notation[m][n]
            stiffness_tensor[j][i][k][l] = voigt_notation[m][n]
            stiffness_tensor[i][j][l][k] = voigt_notation[m][n]
            stiffness_tensor[j][i][l][k] = voigt_notation[m][n]
    return stiffness_tensor


def voigt_notation_to_compliance_tensor(voigt_notation):
    return voigt_notation_to_stiffness_tensor(
        np.divide(voigt_notation, voigt_compliance_factors)
    )


def contract_stresses(stresses):
    """
    Takes a stress tensor in standard (3x3) form
    and returns the Voigt form (6). No factors
    are required to maintain the relationship
    with the corresponding stiffness and strain tensors
    (see contract_strains, which does require
    multiplicative factors).
    """
    return stresses[[0, 1, 2, 1, 0, 0], [0, 1, 2, 2, 2, 1]]


def expand_stresses(stresses):
    """
    Takes a stress tensor in Voigt form (6)
    and returns the standard form (3x3). No factors
    are required to maintain the relationship
    with the corresponding stiffness and strain tensors
    (see contract_strains, which does require
    multiplicative factors).
    """
    return np.array(
        [
            [stresses[0], stresses[5], stresses[4]],
            [stresses[5], stresses[1], stresses[3]],
            [stresses[4], stresses[3], stresses[2]],
        ]
    )


def contract_strains(strains):
    """
    Takes a stress tensor in standard (3x3) form
    and returns the Voigt form (6). Note the factors
    which are required to maintain the relationship
    with the corresponding stiffness and strain tensors.
    """
    # next line creates a copy, not just a view.
    eps = strains[[0, 1, 2, 1, 0, 0], [0, 1, 2, 2, 2, 1]]
    # only overwrites values in eps, not the input array
    eps[3:] *= 2.0
    return eps


def contract_compliances(compliances):
    """
    Takes a compliance tensor in standard (3x3x3x3) form
    and returns the Voigt form (6x6). Note the compliance factors
    which are required to maintain the inverse relationship with the
    corresponding stiffness tensor.
    """
    voigt_notation = np.zeros((6, 6))
    for p in range(6):
        i, j = voigt_index_to_ij(p)
        for q in range(6):
            m, n = voigt_index_to_ij(q)
            voigt_notation[p, q] = compliances[i, j, m, n]
    return np.multiply(voigt_notation, voigt_compliance_factors)


def contract_stiffnesses(stiffnesses):
    """
    Takes a stiffness tensor in standard (3x3x3x3) form
    and returns the Voigt form (6x6).
    """
    voigt_notation = np.zeros((6, 6))
    for p in range(6):
        i, j = voigt_index_to_ij(p)
        for q in range(6):
            m, n = voigt_index_to_ij(q)
            voigt_notation[p, q] = stiffnesses[i, j, m, n]
    return voigt_notation


def voigt_array_from_cijs(cijs, index_lists):
    """
    Takes a list of cijs and a list of list of tuples corresponding to
    the positions of each cij in the Voigt form matrix.
    Note that the indices run from 0--5, not 1--6.
    """
    C = np.zeros([6, 6])
    for i, index_list in enumerate(index_lists):
        for indices in index_list:
            C[indices] = cijs[i]
            C[indices[::-1]] = cijs[i]
    return C
