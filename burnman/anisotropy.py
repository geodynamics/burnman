import numpy as np
import matplotlib.pyplot as plt

def _voigt_index_to_ij(m):
    """
    Returns the ij (or kl) indices of the 
    stiffness tensor which correspond to those 
    of the Voigt notation m (or n).
    """
    if m == 3:
        i=1
        j=2
    elif m == 4:
        i=0
        j=2
    elif m == 5:
        i=0
        j=1
    else:
        i=m
        j=m
    return i, j
            
def voigt_notation_to_stiffness_tensor(voigt_notation):
    """
    Converts a stiffness tensor in Voigt notation 
    to the full fourth rank tensor.

    Parameters
    ----------
    voigt_notation : numpy array
        6 x 6 array

    Returns
    -------
    stiffness_tensor : numpy array
        3 x 3 x 3 x 3 array
    """

    stiffness_tensor = np.zeros([3, 3, 3, 3])
    for m in xrange(6):
        i, j = _voigt_index_to_ij(m)
        for n in xrange(6):
            k, l = _voigt_index_to_ij(n)
            stiffness_tensor[i][j][k][l] = voigt_notation[m][n]
            stiffness_tensor[j][i][k][l] = voigt_notation[m][n]
            stiffness_tensor[i][j][l][k] = voigt_notation[m][n]
            stiffness_tensor[j][i][l][k] = voigt_notation[m][n]
    return stiffness_tensor
    
            
def christoffel_tensor(stiffness_tensor, propagation_direction):
    """
    Computes the Christoffel tensor from an elastic stiffness
    tensor and a propagation direction for a seismic wave.

    # T_ik = C_ijkl n_j n_l

    Parameters
    ----------
    stiffness_tensor : numpy array
        3 x 3 x 3 x 3 array

    Returns
    -------
    Tik : numpy array
        3 x 3 array
    """
    
    Tik = np.tensordot(np.tensordot(stiffness_tensor,
                                    propagation_direction,
                                    axes=([1],[0])),
                       propagation_direction,
                       axes=([2],[0]))
    return Tik

def compliance_tensor(tensor):
    """
    Computes the inverse of a tensor.
    """
    return np.linalg.inv(tensor)


def volume_compressibility(stiffness_tensor):
    """
    Computes the volume compressibility
    from a stiffness tensor

    Parameters
    ----------
    stiffness_tensor : numpy array
        3 x 3 x 3 x 3 array

    Returns
    -------
    beta : float
        volume compressibility
    """
    Sijkl = compliance_tensor(stiffness_tensor)
    Sijkl = voigt_notation_to_stiffness_tensor(Sijkl)
    beta = 0.
    for i in range(3):
        for k in range(3):
            beta += Sijkl[i][i][k][k]
    return beta

def linear_compressibility(stiffness_tensor, direction):
    """
    Computes the linear compressibility
    from a stiffness tensor

    Parameters
    ----------
    stiffness_tensor : numpy array
        3 x 3 x 3 x 3 array
    direction : numpy array
        1D array

    Returns
    -------
    beta : float
        linear compressibility
    """
    Sijkl = compliance_tensor(stiffness_tensor)
    Sijkl = voigt_notation_to_stiffness_tensor(Sijkl)
    
    Sijkk = np.einsum('ijkk', Sijkl)

    beta = np.dot(np.dot(Sijkk, direction),
                        direction)
    return beta

def youngs_modulus(stiffness_tensor, direction):
    Sijkl = compliance_tensor(stiffness_tensor)
    Sijkl = voigt_notation_to_stiffness_tensor(Sijkl)
    S = np.tensordot(np.tensordot(np.tensordot(np.tensordot(Sijkl, direction, axes=([3], [0])),
                             direction, axes=([2], [0])),
                      direction, axes=([1], [0])),
                     direction, axes=([0], [0]))
    return 1./S

def shear_modulus(stiffness_tensor,
                  plane_normal, shear_direction):
    Sijkl = compliance_tensor(stiffness_tensor)
    Sijkl = voigt_notation_to_stiffness_tensor(Sijkl)
    G = np.dot(np.dot(np.dot(np.dot(Sijkl, shear_direction),
                             plane_normal),
                      shear_direction),
               plane_normal)
    return 0.25/G


def poissons_ratio(stiffness_tensor,
                   longitudinal_direction,
                   transverse_direction):
    Sijkl = compliance_tensor(stiffness_tensor)
    Sijkl = voigt_notation_to_stiffness_tensor(Sijkl)
    
    num = np.dot(np.dot(np.dot(np.dot(Sijkl, longitudinal_direction),
                               longitudinal_direction),
                        longitudinal_direction),
                 longitudinal_direction)
    denom = np.dot(np.dot(np.dot(np.dot(Sijkl, transverse_direction),
                                 transverse_direction),
                          transverse_direction),
                   transverse_direction)
    return num/denom

def wave_velocities(stiffness_tensor, propagation_direction, density):
    stiffness_tensor = voigt_notation_to_stiffness_tensor(stiffness_tensor)
    Tik = christoffel_tensor(stiffness_tensor, propagation_direction)

    eigenvalues, eigenvectors = np.linalg.eig(Tik)

    idx = eigenvalues.argsort()[::-1]   
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:,idx]
    
    velocities = np.sqrt(eigenvalues/density)

    return velocities, eigenvectors

