# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import matplotlib.pyplot as plt

import burnman.tools

class SeismicModel:
    """
    Base class for all the seismological models.
    """
    def __init__(self):
        pass
    

        
    def evaluate_all_at(self, depth_list):
        """
        Returns the lists of data for a SeismicModel for the depths provided
            
        Parameters
        ----------
        depth_list : array of floats
            Array of depths (m) to evaluate seismic model at.
            
        Returns
        -------
        pressures : array of floats
            Pressures at given depths in [Pa].
        density :  array of floats
            Densities at given depths in [kg/m^3].
        v_p : array of floats
            P wave velocities at given depths in [m/s].
        v_s : array of floats
            S wave velocities at given depths in [m/s].
        v_phi : array of floats
            Bulk sound velocities at given depths in [m/s].
        """
        pressures = np.array([self.pressure(r) for r in depth_list])
        density = np.array([self.density(r) for r in depth_list])
        v_p = np.array([self.v_p(r) for r in depth_list])
        v_s = np.array([self.v_s(r) for r in depth_list])
        v_phi = np.array([self.v_phi(r) for r in depth_list])
        return pressures, density, v_p, v_s, v_phi

    def _pressure(self, depth):
        raise ValueError, "not implemented"
        return 0

    def _v_p(self, depth):
        raise ValueError, "not implemented"
        return 0

    def _v_s(self, depth):
        raise ValueError, "not implemented"
        return 0
    
    def _v_phi(self, depth):
        raise ValueError, "not implemented"
        return 0

    def _density(self, depth):
        raise ValueError, "not implemented"
        return 0
    
    def _depth(self, pressure):
        raise ValueError, "not implemented"
        return -1


class Model1D(SeismicModel):
    """ 
    This is a base class that gets a 1D seismic model from a table indexed and
    sorted by radius. Fill the tables in the constructor after deriving
    from this class. This class uses :class:`burnman.seismic.SeismicModel`
    
    Format of the input files is = Radius [m], Pressure [Pa], Density [kg/m^3], Vp [m/s], Vs [m/s]
    
    Note: all tables need to be sorted by increasing radius.
    Alternatively, you can also overwrite the _lookup function if you
    want to access with something else.
    """ 
    def __init__(self):
        SeismicModel.__init__(self)
        self.table_radius = []
        self.table_pressure = []
        self.table_density = []
        self.table_vp = []
        self.table_vs = []
        self.earth_radius = 6371.0e3
        
    def internal_depth_list(self):
        '''
        Returns
        -------
        depths : array of floats
            Depths [m].
        '''
        
        return (self.earth_radius - self.table_radius)[::-1] #radius is sorted in increasing order, so we need to reverse the depth list

    def pressure(self, depth):
        '''
        Parameters
        ----------
        depth_list : array of floats
            Array of depths (m) to evaluate seismic model at.
            
        Returns
        -------
        pressures : array of floats
            Pressures at given depths in [Pa].
        '''
        return self._lookup(depth, self.table_pressure)

    def v_p(self, depth):
        '''
        Parameters
        ----------
        depth_list : array of floats
            Array of depths (m) to evaluate seismic model at.
            
        Returns
        -------
        v_p : array of floats
            P wave velocities at given depths in [m/s].
        '''
        return self._lookup(depth, self.table_vp)

    def v_s(self, depth):
        '''
        Parameters
        ----------
        depth_list : array of floats
            Array of depths (m) to evaluate seismic model at.
            
        Returns
        -------
        v_s : array of floats
            S wave velocities at given depths in [m/s].
        '''
        return self._lookup(depth, self.table_vs)

    def v_phi(self, depth):
        '''
        Parameters
        ----------
        depth_list : array of floats
            Array of depths (m) to evaluate seismic model at.
            
        Returns
        -------
        v_phi : array of floats
            bulk sound wave velocities at given depths in [m/s].
            '''
        
        v_s=self.v_s(depth)
        v_p=self.v_p(depth)
        return np.sqrt(v_p*v_p-4./3.*v_s*v_s)

    def density(self, depth):
        '''
        Parameters
        ----------
        depth_list : array of floats
        Array of depths (m) to evaluate seismic model at.
            
        Returns
        -------
        density : array of floats
            densities at given depths in [kg/m^3].
        '''
        return self._lookup(depth, self.table_density)        

    def depth(self, pressure):
        '''
        Parameters
        ----------
        pressures : array of floats
            Pressures at given depths in [Pa].
            
        Returns
        -------
        depth_list : array of floats
            Array of depths (m) to evaluate seismic model at.
        '''
        radius = burnman.tools.lookup_and_interpolate(self.table_pressure[::-1], self.table_radius[::-1], pressure)
        return self.earth_radius - radius

    def _lookup(self, depth, value_table):
        radius = self.earth_radius - depth
        return burnman.tools.lookup_and_interpolate(self.table_radius, value_table, radius)    
    

class PREM(Model1D):
    """ 
    Reads  PREM (1s) (input_seismic/prem_table.txt, Dziewonski & Anderson 1981). 
    See also :class:`burnman.seismic.Model1D`.
    """
    def __init__(self):
        Model1D.__init__(self)
        table = burnman.tools.read_table("input_seismic/prem_table.txt") # radius, pressure, density, v_p, v_s
        table = np.array(table)
        self.table_radius = table[:,0]
        self.table_pressure = table[:,1] 
        self.table_density = table[:,2]
        self.table_vp = table[:,3]
        self.table_vs = table[:,4]

    def grav(self,depths):
        '''
        Parameters
        ----------
        depth_list : array of floats
            Array of depths (m) to evaluate seismic model at.
            
        Returns
        -------
        grav : array of floats
            Gravity for PREM model for given depths in [m/s^2]
        '''
        
        table = burnman.tools.read_table("input_seismic/grav_for_PREM.txt") # radius, g
        table = np.array(table)
        table_rad = table[:,0]
        table_g = table[:,1]
        return np.interp(self.earth_radius-depths, table_rad,table_g)


class Slow(Model1D):
    """ 
    Inserts the mean profiles for slower regions in the lower mantle (Lekic et al. 2012). 
    We stitch together tables 'input_seismic/prem_lowermantle.txt', 'input_seismic/swave_slow.txt', 'input_seismic/pwave_slow.txt').
    See also :class:`burnman.seismic.Model1D`.
    """
    def __init__(self):
        Model1D.__init__(self)

        table = burnman.tools.read_table("input_seismic/prem_lowermantle.txt")#data is: radius pressure density V_p V_s Q_K Q_G
        table = np.array(table)
        table[:,0] = table[:,0]  
        table2 = burnman.tools.read_table("input_seismic/swave_slow.txt")
        table2 = np.array(table2)
        table3 = burnman.tools.read_table("input_seismic/pwave_slow.txt")
        table3 = np.array(table3)

        min_radius = self.earth_radius-max(table2[:,0])
        max_radius = self.earth_radius-min(table2[:,0])
        
        table=np.array(filter(lambda x: (x[0]>=min_radius and x[0]<=max_radius), table))

        assert(len(table) == len(table2))
        assert(len(table) == len(table3))

        self.table_radius = table[:,0]
        self.table_pressure = table[:,1] 
        self.table_density = table[:,2]
        self.table_vp = table3[:,1] 
        self.table_vs = table2[:,1] 
       

class Fast(Model1D):
    """ 
    Inserts the mean profiles for faster regions in the lower mantle (Lekic et al. 2012). 
    We stitch together tables 'input_seismic/prem_lowermantle.txt', 'input_seismic/swave_fast.txt', 'input_seismic/pwave_fast.txt').
    See also :class:`burnman.seismic.Model1D`.
    """
    def __init__(self):
        Model1D.__init__(self)

        table = burnman.tools.read_table("input_seismic/prem_lowermantle.txt")#data is: radius pressure density V_p V_s Q_K Q_G
        table = np.array(table)
        table[:,0] = table[:,0]  
        table2 = burnman.tools.read_table("input_seismic/swave_fast.txt")
        table2 = np.array(table2)
        table3 = burnman.tools.read_table("input_seismic/pwave_fast.txt")
        table3 = np.array(table3)

        min_radius = self.earth_radius-max(table2[:,0])
        max_radius = self.earth_radius-min(table2[:,0])
        
        table=np.array(filter(lambda x: (x[0]>=min_radius and x[0]<=max_radius), table))

        assert(len(table) == len(table2))
        assert(len(table) == len(table3))

        self.table_radius = table[:,0]
        self.table_pressure = table[:,1] 
        self.table_density = table[:,2]
        self.table_vp = table3[:,1] 
        self.table_vs = table2[:,1] 




def attenuation_correction(v_p,v_s,v_phi,Qs,Qphi):
    """ 
    Applies the attenuation correction following Matas et al. (2007), page 4. This is simplified, and there is also currently no 1D Q model implemented. The correction, however, only slightly reduces the velocities, and can be ignored for our current applications. Arguably, it might not be as relevant when comparing computations to PREM for periods of 1s as is implemented here.
    Called from :func:`burnman.main.apply_attenuation_correction`
    
    Parameters
    ----------
    v_p : float
        P wave velocity in [m/s].
    v_s : float
        S wave velocitiy in [m/s].
    v_phi : float
        Bulk sound velocity in [m/s].
    Qs : float
        shear quality factor [dimensionless]
    Qphi: float
        bulk quality factor [dimensionless]
        
    
    Returns
    -------
    v_p : float
    corrected P wave velocity in [m/s].
    v_s : float
    corrected S wave velocitiy in [m/s].
    v_phi : float
    corrected Bulk sound velocity in [m/s].



    """
    beta = 0.3 # Matas et al. (2007) page 4
    Qp  = 3./4.*pow((v_p/v_s),2.)*Qs    # Matas et al. (2007) page 4


    cot=1./np.tan(beta*np.pi/2.)
    v_p  = v_p*(1.-1./2.*cot*1./Qp)    # Matas et al. (2007) page 1
    v_s  = v_s*(1.-1./2.*cot*1./Qs)
    v_phi= v_phi*(1.-1./2.*cot*1./Qphi)
    return v_p, v_s, v_phi
    
# shared variable of prem, so that other routines do not need to create
# prem over and over. See geotherm for example.
prem_model = PREM()

if __name__ == "__main__":
    #create a seismic dataset from prem:
    s=PREM()
    depths = s.internal_depth_list()
    pressures, density, v_p, v_s, v_phi = s.evaluate_all_at(depths)

    # plot vs and vp and v_phi (note that v_phi is computed!)
    plt.plot(depths/1.e3,v_p/1.e3,'+-r', label='v_p')
    plt.plot(depths/1.e3,v_s/1.e3,'+-b', label='v_s')
    plt.plot(depths/1.e3,v_phi/1.e3,'--g', label='v_phi')
    plt.legend()
    plt.xlabel('depth in km')
    plt.ylabel('km/s')
    plt.show()

    s1=PREM()
    depths=s1.internal_depth_list()
    pressures, density, v_p, v_s, v_phi = s1.evaluate_all_at(depths)
    plt.plot(depths/1.e3,v_p/1.e3,'+-r', label='v_p')
    plt.plot(depths/1.e3,v_s/1.e3,'+-b', label='v_s')
    plt.plot(depths/1.e3,v_phi/1.e3,'--g', label='v_phi')


    
