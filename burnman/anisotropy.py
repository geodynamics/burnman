# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt

from .tools import normalize

class AnisotropicMaterial(object):
    """
    A class that represents an anisotropic elastic material. This class 
    is initialised with a set of elastic constants and a density. It can
    then be interrogated to find the values of different properties, 
    such as bounds on seismic velocities. There are also several functions
    which can be called to calculate properties along directions oriented
    with respect to the elastic tensor.

    Initialization is with a density and either 
    a) a set of independent elastic constants and a crystal system, or
    b) a full stiffness tensor in Voigt notation

    If initialization is via option (a), the number and order 
    of the constants is dependent on the crystal system:
    'isotropic': C12, C44 (i.e. lambda and mu, the Lame parameters)
    'cubic': C11, C12, C44
    'hexagonal': C11, C12, C13, C33, C44
    'tetragonal I': C11, C12, C13, C33, C44, C66
    'tetragonal II': C11, C12, C13, C16, C33, C44, C66
    'rhombohedral I': C11, C12, C13, C14, C33, C44, C66
    'rhombohedral II': C11, C12, C13, C14, C15, C33, C44, C66
    'orthorhombic': C11, C12, C13, C22, C23, C33, C44, C55, C66
    'monoclinic': C11, C12, C13, C15, C22, C23, 
                  C25, C33, C35, C44, C46, C55, C66
    'triclinic': Cij, where 1<=i<=6 and i<=j<=6 
    
    See :cite:`Mainprice2011` Geological Society of London Special Publication 
    and https://materialsproject.org/wiki/index.php/Elasticity_calculations
    for mathematical descriptions of each function.
    """
        
    def __init__(self, cijs, rho, crystal_system=None):
        if crystal_system is not None:
            self.stiffness_tensor = self._cijs_to_voigt(cijs, crystal_system)
        else:
            self.stiffness_tensor = np.array(cijs)
        
        assert self.stiffness_tensor.shape == (6, 6), 'stiffness_tensor must be in Voigt notation (6x6)'
        assert np.allclose(self.stiffness_tensor.T,
                           self.stiffness_tensor), 'stiffness_tensor must be symmetric'
        
        self.full_stiffness_tensor = self._voigt_notation_to_stiffness_tensor(self.stiffness_tensor)

        self.compliance_tensor = np.linalg.inv(self.stiffness_tensor)
        
        block = np.array(np.bmat( [[[[1.]*3]*3, [[2.]*3]*3], [[[2.]*3]*3, [[4.]*3]*3]] ))
        self.full_compliance_tensor = self._voigt_notation_to_stiffness_tensor(np.divide(self.compliance_tensor, block))
        self.rho = rho

    def _cijs_to_voigt(self, cijs, crystal_system):
        """
        Converts individual elastic tensors (cijs) to
        the stiffness tensor in Voigt notation (6x6)
        based on the crystal system. A list of crystal systems 
        is provided in the base class description.
        """
        
        if crystal_system=='isotropic':
            assert len(cijs) == 2
            cijs = list(cijs)
            cijs.insert(0, cijs[0] + 2.*cijs[1]) # C11 = C12 + 2C44
            index_lists = [[(0, 0), (1, 1), (2, 2)], # C11
                           [(0, 1), (0, 2), (1, 2)], # C12
                           [(3, 3), (4, 4), (5, 5)]] # C44
            
        elif crystal_system=='cubic':
            assert len(cijs) == 3
            index_lists = [[(0, 0), (1, 1), (2, 2)], # C11
                           [(0, 1), (0, 2), (1, 2)], # C12
                           [(3, 3), (4, 4), (5, 5)]] # C44
            
        elif crystal_system=='hexagonal' or crystal_system=='tetragonal I':
            if len(cijs) == 5: #for hexagonal, C66 = (C11-C12)/2.
                cijs = list(cijs)
                cijs.append((cijs[0] - cijs[1])/2.)
            # tetragonal_i corresponds to Laue class 4/mmm
            assert len(cijs) == 6
            index_lists = [[(0, 0), (1, 1)], # C11
                           [(0, 1)], # C12
                           [(0, 2), (1, 2)], # C13
                           [(2, 2)], # C33
                           [(3, 3), (4, 4)], # C44
                           [(5, 5)]] # C66
            
        elif crystal_system=='tetragonal II':
            # tetragonal_ii corresponds to Laue class 4/m
            assert len(cijs) == 7
            cijs = list(cijs)
            cijs.insert(4, cij[3]) # C26 = -C16
            index_lists = [[(0, 0), (1, 1)], # C11
                           [(0, 1)], # C12
                           [(0, 2), (1, 2)], # C13
                           [(0, 5)], # C16
                           [(1, 5)], # C26
                           [(2, 2)], # C33
                           [(3, 3), (4, 4)], # C44
                           [(5, 5)]] # C66
            
        elif crystal_system=='rhombohedral I':
            # rhombohedral_i corresponds to Laue class \bar{3}m
            assert len(cijs) == 7
            cijs = list(cijs)
            cijs.insert(4, cij[3]) # C24 = -C14
            index_lists = [[(0, 0), (1, 1)], # C11
                           [(0, 1)], # C12
                           [(0, 2), (1, 2)], # C13
                           [(0, 3), (4, 5)], # C14
                           [(1, 3)], # C24
                           [(2, 2)], # C33
                           [(3, 3), (4, 4)], # C44
                           [(5, 5)]] # C66
            
        elif crystal_system=='rhombohedral II':
            # rhombohedral_ii corresponds to Laue class \bar{3}
            assert len(cijs) == 8
            cijs = list(cijs)
            cijs.insert(4, cij[3]) # C24 = -C14
            cijs.insert(6, cij[5]) # C25 = -C15
            index_lists = [[(0, 0), (1, 1)], # C11
                           [(0, 1)], # C12
                           [(0, 2), (1, 2)], # C13
                           [(0, 3), (4, 5)], # C14
                           [(1, 3)], # C24
                           [(0, 4)], # C15
                           [(1, 4), (3, 5)], # C25
                           [(2, 2)], # C33
                           [(3, 3), (4, 4)], # C44
                           [(5, 5)]] # C66
            
        elif crystal_system=='orthorhombic':
            assert len(cijs) == 9
            index_lists = [[(0, 0)], # C11
                           [(0, 1)], # C12
                           [(0, 2)], # C13
                           [(1, 1)], # C22
                           [(1, 2)], # C23
                           [(2, 2)], # C33
                           [(3, 3)], # C44
                           [(4, 4)], # C55
                           [(5, 5)]] # C66
            
        elif crystal_system=='monoclinic':
            assert len(cijs) == 13
            index_lists = [[(0, 0)], # C11
                           [(0, 1)], # C12
                           [(0, 2)], # C13
                           [(0, 4)], # C15
                           [(1, 1)], # C22
                           [(1, 2)], # C23
                           [(1, 4)], # C25
                           [(2, 2)], # C33
                           [(2, 4)], # C35
                           [(3, 3)], # C44
                           [(3, 5)], # C46
                           [(4, 4)], # C55
                           [(5, 5)]] # C66
            
        elif crystal_system=='triclinic':
            assert len(cijs) == 21
            index_lists=[[(i, j)] for i in range(6) for j in range(i, 6)]
            
        else:
            raise Exception('Crystal system not recognised. Must be one of: '
                            'isotropic, cubic, hexagonal, tetragonal I, '
                            'tetragonal II, rhombohedral I, rhombohedral II, '
                            'orthorhombic, monoclinic or triclinic.')    
            
        
        C = np.zeros([6, 6])
        for i, index_list in enumerate(index_lists):
            for indices in index_list:
                C[indices] = cijs[i]
                C[indices[::-1]] = cijs[i]
        return C
        
    def _voigt_index_to_ij(self, m):
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

    def _voigt_notation_to_stiffness_tensor(self, voigt_notation):
        """
        Converts a stiffness tensor in Voigt notation (6x6 matrix)
        to the full fourth rank tensor (3x3x3x3 matrix).
        """
        stiffness_tensor = np.zeros([3, 3, 3, 3])
        for m in range(6):
            i, j = self._voigt_index_to_ij(m)
            for n in range(6):
                k, l = self._voigt_index_to_ij(n)
                stiffness_tensor[i][j][k][l] = voigt_notation[m][n]
                stiffness_tensor[j][i][k][l] = voigt_notation[m][n]
                stiffness_tensor[i][j][l][k] = voigt_notation[m][n]
                stiffness_tensor[j][i][l][k] = voigt_notation[m][n]
        return stiffness_tensor
    
    @property
    def bulk_modulus_voigt(self):
        """
        Computes the bulk modulus (Voigt bound)
        """
        K = np.sum([[self.stiffness_tensor[i][k] for k in range(3)] for i in range(3)])/9.
        return K
    
    @property
    def bulk_modulus_reuss(self):
        """
        Computes the bulk modulus (Reuss bound)
        """
        beta = np.sum([[self.compliance_tensor[i][k] for k in range(3)] for i in range(3)])
        return 1./beta

    @property
    def bulk_modulus_vrh(self):
        """
        Computes the bulk modulus (Voigt-Reuss-Hill average)
        """
        return 0.5*(self.bulk_modulus_voigt + self.bulk_modulus_reuss)
    
    @property
    def shear_modulus_voigt(self):
        """
        Computes the shear modulus (Voigt bound)
        """
        G = ( np.sum([self.stiffness_tensor[i][i] for i in [0, 1, 2]]) +
              np.sum([self.stiffness_tensor[i][i] for i in [3, 4, 5]])*3. -
              ( self.stiffness_tensor[0][1] +
                self.stiffness_tensor[1][2] +
                self.stiffness_tensor[2][0] )) / 15.
        return G
    
    @property
    def shear_modulus_reuss(self):
        """
        Computes the shear modulus (Reuss bound)
        """
        beta =  ( np.sum([self.compliance_tensor[i][i] for i in [0, 1, 2]])*4. +
                  np.sum([self.compliance_tensor[i][i] for i in [3, 4, 5]])*3. -
                  ( self.compliance_tensor[0][1] +
                    self.compliance_tensor[1][2] + 
                    self.compliance_tensor[2][0])*4. ) / 15.
        return 1./beta
    
    @property
    def shear_modulus_vrh(self):
        """
        Computes the shear modulus (Voigt-Reuss-Hill average)
        """
        return 0.5*(self.shear_modulus_voigt + self.shear_modulus_reuss)

    @property
    def universal_elastic_anisotropy(self):
        """
        Compute the universal elastic anisotropy
        """
        return ( 5.*(self.shear_modulus_voigt/self.shear_modulus_reuss) +
                 (self.bulk_modulus_voigt/self.bulk_modulus_reuss) - 6. )

    @property
    def isotropic_poisson_ratio(self):
        """
        Compute mu, the isotropic Poisson ratio
        (a description of the laterial response to loading)
        """
        return ( (3.*self.bulk_modulus_vrh - 2.*self.shear_modulus_vrh) /
                 (6.*self.bulk_modulus_vrh + 2.*self.shear_modulus_vrh) )

    def christoffel_tensor(self, propagation_direction):
        """
        Computes the Christoffel tensor from an elastic stiffness
        tensor and a propagation direction for a seismic wave
        relative to the stiffness tensor

        T_ik = C_ijkl n_j n_l
        """
        propagation_direction = normalize(propagation_direction)
        Tik = np.tensordot(np.tensordot(self.full_stiffness_tensor,
                                        propagation_direction,
                                        axes=([1],[0])),
                           propagation_direction,
                           axes=([2],[0]))
        return Tik

    def linear_compressibility(self, direction):
        """
        Computes the linear compressibility in a given direction 
        relative to the stiffness tensor
        """
        direction = normalize(direction)
        Sijkk = np.einsum('ijkk', self.full_compliance_tensor)
        beta = Sijkk.dot(direction).dot(direction)
        return beta
    
    def youngs_modulus(self, direction):
        """
        Computes the Youngs modulus in a given direction 
        relative to the stiffness tensor
        """
        direction = normalize(direction)
        Sijkl = self.full_compliance_tensor
        S = Sijkl.dot(direction).dot(direction).dot(direction).dot(direction)
        return 1./S

    def shear_modulus(self, plane_normal, shear_direction):
        """
        Computes the shear modulus on a plane in a given 
        shear direction relative to the stiffness tensor
        """
        plane_normal = normalize(plane_normal)
        shear_direction = normalize(shear_direction)
        
        assert np.abs(plane_normal.dot(shear_direction)) < np.finfo(np.float).eps, 'plane_normal and shear_direction must be orthogonal'
        Sijkl = self.full_compliance_tensor
        G = Sijkl.dot(shear_direction).dot(plane_normal).dot(shear_direction).dot(plane_normal)
        return 0.25/G
    
    def poissons_ratio(self,
                       axial_direction,
                       lateral_direction):
        """
        Computes the poisson ratio given loading and response 
        directions relative to the stiffness tensor
        """
        
        axial_direction = normalize(axial_direction)
        lateral_direction = normalize(lateral_direction)
        assert np.abs(axial_direction.dot(lateral_direction)) < np.finfo(np.float).eps, 'axial_direction and lateral_direction must be orthogonal'
        
        Sijkl = self.full_compliance_tensor
        x = axial_direction
        y = lateral_direction
        nu = -(Sijkl.dot(y).dot(y).dot(x).dot(x) / 
               Sijkl.dot(x).dot(x).dot(x).dot(x) )
        return nu
    
    def wave_velocities(self, propagation_direction):
        """
        Computes the compressional wave velocity, and two 
        shear wave velocities in a given propagation direction

        Returns two lists, containing the wave speeds and 
        directions of particle motion relative to the stiffness tensor
        """
        propagation_direction = normalize(propagation_direction)
        
        Tik = self.christoffel_tensor(propagation_direction)

        eigenvalues, eigenvectors = np.linalg.eig(Tik)

        idx = eigenvalues.argsort()[::-1]   
        eigenvalues = np.real(eigenvalues[idx])        
        eigenvectors = eigenvectors[:,idx]
        velocities = np.sqrt(eigenvalues/self.rho)

        return velocities, eigenvectors

    def plot_velocities(self):
        """
        Makes colour plots of:
        Compressional wave velocity: Vp
        Anisotropy: (Vs1 - Vs2)/(Vs1 + Vs2)
        Vp/Vs1
        linear compressibility: beta
        Youngs Modulus: E
        """

        try:
            plt.style.use('ggplot')
            plt.rcParams['axes.facecolor'] = 'white'
            plt.rcParams['axes.edgecolor'] = 'black'
            plt.rcParams['figure.figsize'] = 16, 10 # inches
        except:
            pass
        
        zeniths = np.linspace(np.pi/2., np.pi, 31)
        azimuths = np.linspace(0., 2.*np.pi, 91)
        Rs = np.sin(zeniths)/(1. - np.cos(zeniths))
        r, theta = np.meshgrid(Rs, azimuths)
        
        vps = np.empty_like(r)
        vs1s = np.empty_like(r)
        vs2s = np.empty_like(r)
        betas = np.empty_like(r)
        Es = np.empty_like(r)
        for i, az in enumerate(azimuths):
            for j, phi in enumerate(zeniths):
                d = np.array([np.cos(az)*np.sin(phi), np.sin(az)*np.sin(phi), -np.cos(phi)]) # change_hemispheres
                velocities = self.wave_velocities(d)
                betas[i][j] = self.linear_compressibility(d)
                Es[i][j] = self.youngs_modulus(d)
                vps[i][j] = velocities[0][0]
                vs1s[i][j] = velocities[0][1]
                vs2s[i][j] = velocities[0][2]
                
        fig = plt.figure()
        names = ['Vp (km/s)', 'Vs1 (km/s)', 'Vp/Vs1', 'S-wave anisotropy (%)', 'linear beta (GPa^-1)', 'Youngs Modulus (GPa)']
        items = [vps/1000., vs1s/1000., vps/vs1s, 200.*(vs1s - vs2s)/(vs1s + vs2s), betas*1.e9, Es/1.e9]
        ax = []
        im = []
        ndivs = 100
        for i, item in enumerate(items):
            ax.append(fig.add_subplot(2, 3, i+1, projection='polar'))
            ax[i].set_yticks([100])
            ax[i].set_title(names[i])

            vmin = np.min(item)
            vmax = np.max(item)
            spacing = np.power(10., np.floor(np.log10(vmax - vmin)))
            nt = int((vmax - vmin - vmax%spacing + vmin%spacing)/spacing)
            if nt == 1:
                spacing = spacing/4.
            elif nt < 4:
                spacing = spacing/2.
            elif nt > 8:
                spacing = spacing*2.
                
            tmin = vmin + (spacing - vmin%spacing)
            tmax = vmax - vmax%spacing
            nt = int((tmax - tmin)/spacing + 1)
            
            ticks = np.linspace(tmin, tmax, nt)
            im.append(ax[i].contourf(theta, r, item, ndivs, cmap=plt.cm.jet_r, vmin=vmin, vmax=vmax))
            lines = ax[i].contour(theta, r, item, ticks, colors=('black',), linewidths=(1,))
            
            cbar = fig.colorbar(im[i], ax=ax[i], ticks=ticks)
            cbar.add_lines(lines)

        plt.tight_layout()
        plt.show()
