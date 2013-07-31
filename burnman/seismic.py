# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import tools
import matplotlib.pyplot as plt
import math

class seismic_data:
    """
    base class for all seismic models
    """
    def __init__(self):
        1 # how else to make an empty function?
        
    # depth in km
    def internal_depth_list(self):
        """ returns a sorted list of depths where this seismic data is computed at """
        return np.arange(0.,6000.0, 100.0)
        
    def evaluate_all_at(self, depth_list):
	""" returns pressure[Pa], density[kg/m^3], Vp[m/s], Vs[m/s] and Vphi[m/s] for a list of depths[m] """
        pressures = np.array([self.pressure(r) for r in depth_list])
        density = np.array([self.density(r) for r in depth_list])
        v_p = np.array([self.v_p(r) for r in depth_list])
        v_s = np.array([self.v_s(r) for r in depth_list])
        v_phi = np.array([self.v_phi(r) for r in depth_list])
        return pressures, density, v_p, v_s, v_phi

    def pressure(self, depth):
        raise ValueError, "not implemented"
        return 0

    def v_p(self, depth):
        raise ValueError, "not implemented"
        return 0

    def v_s(self, depth):
        raise ValueError, "not implemented"
        return 0
    
    def v_phi(self, depth):
        raise ValueError, "not implemented"
        return 0

    def density(self, depth):
        raise ValueError, "not implemented"
        return 0
    
    def depth(self, pressure):
        raise ValueError, "not implemented"
        return -1


class radiustable(seismic_data):
    """ 
    this is a base class that gets the information from a table indexed and
    sorted by radius. Fill the tables in the constructor after deriving
    from this class.
    Note: all tables need to be sorted by increasing radius.
    Alternatively, you can also overwrite the _lookup function if you
    want to access with something else like depth than radius.
    """ 
    def __init__(self):
        seismic_data.__init__(self)
        self.table_radius = []
        self.table_pressure = []
        self.table_density = []
        self.table_vp = []
        self.table_vs = []
        self.earth_radius = 6371.0e3
        
    def internal_depth_list(self):
        return (self.earth_radius - self.table_radius)[::-1] #radius is sorted in increasing order, so we need to reverse the depth list

    def pressure(self, depth):
        return self._lookup(depth, self.table_pressure)

    def v_p(self, depth):
        return self._lookup(depth, self.table_vp)

    def v_s(self, depth):
        return self._lookup(depth, self.table_vs)

    def v_phi(self, depth):
        v_s=self.v_s(depth)
        v_p=self.v_p(depth)
        return math.sqrt(v_p*v_p-4./3.*v_s*v_s)

    def density(self, depth):
        return self._lookup(depth, self.table_density)        

    def depth(self, pressure):
        radius = tools.lookup_and_interpolate(self.table_pressure[::-1], self.table_radius[::-1], pressure)
        return self.earth_radius - radius

    def _lookup(self, depth, value_table):
        radius = self.earth_radius - depth
        return tools.lookup_and_interpolate(self.table_radius, value_table, radius)    
    

class prem(radiustable):
    """ 
    reads in the table for PREM (input_seismic/prem_table.txt) using the base class radiustable
    """
    def __init__(self):
        radiustable.__init__(self)
        table = tools.read_table("input_seismic/prem_table.txt") # radius, pressure, density, v_p, v_s
        table = np.array(table)
        self.table_radius = table[:,0]
        self.table_pressure = table[:,1] 
        self.table_density = table[:,2]
        self.table_vp = table[:,3]
        self.table_vs = table[:,4]


class slow(radiustable):
    """ 
    Inserts the mean profiles for slower regions in the lower mantle (Lekic et al. 2012). 
    We need to stitch together three tables. Note that prem_lowermantle has a wider range, 
    so we cut away rows at the top and bottom. Interpolation is not necessary, 
    because all tables where generated with at the same depths 
    """
    def __init__(self):
        radiustable.__init__(self)

        table = tools.read_table("input_seismic/prem_lowermantle.txt")#data is: radius pressure density V_p V_s Q_K Q_G
        table = np.array(table)
        table[:,0] = table[:,0]  
        table2 = tools.read_table("input_seismic/swave_slow.txt")
        table2 = np.array(table2)
        table3 = tools.read_table("input_seismic/pwave_slow.txt")
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
       

class fast(radiustable):
    """ 
    Inserts the mean profiles for faster regions in the lower mantle (Lekic et al. 2012). 
    We need to stitch together three tables. Note that prem_lowermantle has a wider range, 
    so we cut away rows at the top and bottom. Interpolation is not necessary, 
    because all tables where generated with at the same depths 
    """
    def __init__(self):
        radiustable.__init__(self)

        table = tools.read_table("input_seismic/prem_lowermantle.txt")#data is: radius pressure density V_p V_s Q_K Q_G
        table = np.array(table)
        table[:,0] = table[:,0]  
        table2 = tools.read_table("input_seismic/swave_fast.txt")
        table2 = np.array(table2)
        table3 = tools.read_table("input_seismic/pwave_fast.txt")
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



# this uses prem_lowermantle table
class prem_test(radiustable):
    def __init__(self):
        radiustable.__init__(self)

        table = tools.read_table("input_seismic/prem_lowermantle.txt")#data is: radius pressure density V_p V_s Q_K Q_G
        table = np.array(table)
        self.table_radius = table[:,0]   
        self.table_pressure = table[:,1] 
        self.table_density = table[:,2]
        self.table_vp = table[:,3] 
        self.table_vs = table[:,4] 

def attenuation_correction(v_p,v_s,v_phi,Qs,Qphi):
    """ 
    Applies the attenuation correction following Matas et al. (2007), page 4. This is a minor effect on the velocities 
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
prem_model = prem()

if __name__ == "__main__":
    #create a seismic dataset from prem:
    s=prem()
    depths = s.internal_depth_list()
    pressures, density, v_p, v_s, v_phi = s.evaluate_all_at(depths)
    print depths, pressures, density, v_p, v_s, v_phi


    # specify where we want to evaluate, here we map from pressure to depth, because we can
    #p = np.arange(1e9,360e9,5e9)
    #depths = map(s.depth, p) 
    #we could also just specify some depth levels directly like this:
    #depths = np.arange(35e3,5600e3,100e3)

    #now evaluate everything at the given depths levels (using interpolation)
    #pressures, density, v_p, v_s, v_phi = s.evaluate_all_at(depths)

    # plot vs and vp and v_phi (note that v_phi is computed!)
    plt.plot(depths/1.e3,v_p/1.e3,'+-r', label='v_p')
    plt.plot(depths/1.e3,v_s/1.e3,'+-b', label='v_s')
    plt.plot(depths/1.e3,v_phi/1.e3,'--g', label='v_phi')
    plt.legend()
    plt.xlabel('depth in km')
    plt.ylabel('km/s')
    plt.show()

    s1=prem()
    depths=s1.internal_depth_list()
    pressures, density, v_p, v_s, v_phi = s1.evaluate_all_at(depths)
    plt.plot(depths/1.e3,v_p/1.e3,'+-r', label='v_p')
    plt.plot(depths/1.e3,v_s/1.e3,'+-b', label='v_s')
    plt.plot(depths/1.e3,v_phi/1.e3,'--g', label='v_phi')

    s2=prem_test()
    depths=s2.internal_depth_list()
    pressures, density, v_p, v_s, v_phi = s2.evaluate_all_at(depths)
    plt.plot(depths,v_p/1.e3,'x-r', label='v_p')
    plt.plot(depths,v_s/1.e3,'x-b', label='v_s')
    plt.plot(depths,v_phi/1.e3,'x-g', label='v_phi')
    plt.show()
    
