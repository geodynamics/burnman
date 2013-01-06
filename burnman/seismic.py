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
        pressures = np.array([self.pressure(r) for r in depth_list])
        density = np.array([self.density(r) for r in depth_list])
        v_p = np.array([self.v_p(r) for r in depth_list])
        v_s = np.array([self.v_s(r) for r in depth_list])
        v_phi = np.array([self.v_phi(r) for r in depth_list])
        return pressures, density, v_p, v_s, v_phi

    # in Pa, depth in km
    def pressure(self, depth):
        raise ValueError, "not implemented"
        return 0

    #in km/s, depth in km
    def v_p(self, depth):
        raise ValueError, "not implemented"
        return 0

    #in km/s, depth in km
    def v_s(self, depth):
        raise ValueError, "not implemented"
        return 0
    
    #in km/s, depth in km
    def v_phi(self, depth):
        raise ValueError, "not implemented"
        return 0

    #in g/cc, depth in km
    def density(self, depth):
        raise ValueError, "not implemented"
        return 0
    
    #look up depth from given pressure in Pa
    def depth(self, pressure):
        raise ValueError, "not implemented"
        return -1

    #constants for attenuation correction
    #returns (Qk, Qmu) with Qk: factor for v_phi, Qmu: factor for v_s
    def attenuation(self, pressure):
        return 10000., 355.

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
        self.earth_radius = 6371.0
        
    def internal_depth_list(self):
        return (self.earth_radius - self.table_radius)[::-1] #radius is sorted in increasing order, so we need to reverse the depth list

    def pressure(self, depth):
        return self._lookup(depth, self.table_pressure)

    def v_p(self, depth):
        return self._lookup(depth, self.table_vp)

    def v_s(self, depth):
        return self._lookup(depth, self.table_vs)

    def v_phi(self, depth):
        vs=self.v_s(depth)
        vp=self.v_p(depth)
        return math.sqrt(vp*vp-4./3.*vs*vs)

    def density(self, depth):
        return self._lookup(depth, self.table_density)        

    def depth(self, pressure):
        radius = tools.lookup_and_interpolate(self.table_pressure[::-1], self.table_radius[::-1], pressure)
	return self.earth_radius - radius

    def _lookup(self, depth, value_table):
        radius = self.earth_radius - depth
        return tools.lookup_and_interpolate(self.table_radius, value_table, radius)	
    

class prem(radiustable):
    def __init__(self):
        radiustable.__init__(self)
        table = tools.read_table("data/prem_table.txt") # radius, pressure, density, vp, vs
        table = np.array(table)
        self.table_radius = table[:,0]
        self.table_pressure = table[:,1] * 0.1 * 1e9 # convert kbar to GPa to Pa
        self.table_density = table[:,2]
        self.table_vp = table[:,3]
        self.table_vs = table[:,4]


class slow(radiustable):
    def __init__(self):
        radiustable.__init__(self)

        #we need to stitch together three tables. Note that prem_withQ has a wider range, so we cut
        #away rows at the top and bottom. Interpolation is not necessary, because all tables
        #where generated with at the same depths

        table = tools.read_table("data/prem_withQ.txt")#data is: radius pressure density V_p V_s Q_K Q_mu
        table = np.array(table)
        table[:,0] = table[:,0] / 1e3 #convert to km 
        table2 = tools.read_table("data/swave_slow.txt")
        table2 = np.array(table2)
        table3 = tools.read_table("data/pwave_slow.txt")
        table3 = np.array(table3)

        min_radius = self.earth_radius-max(table2[:,0])/1e3
        max_radius = self.earth_radius-min(table2[:,0])/1e3
        
        table=np.array(filter(lambda x: (x[0]>=min_radius and x[0]<=max_radius), table))

        assert(len(table) == len(table2))
        assert(len(table) == len(table3))

        self.table_radius = table[:,0]
        self.table_pressure = table[:,1] * 1e9
        self.table_density = table[:,2]
        self.table_vp = table3[:,1] / 1e3
        self.table_vs = table2[:,1] / 1e3
       

class fast(radiustable):
    def __init__(self):
        radiustable.__init__(self)

        #we need to stitch together three tables. Note that prem_withQ has a wider range, so we cut
        #away rows at the top and bottom. Interpolation is not necessary, because all tables
        #where generated with at the same depths

        table = tools.read_table("data/prem_withQ.txt")#data is: radius pressure density V_p V_s Q_K Q_mu
        table = np.array(table)
        table[:,0] = table[:,0] / 1e3 #convert to km 
        table2 = tools.read_table("data/swave_fast.txt")
        table2 = np.array(table2)
        table3 = tools.read_table("data/pwave_fast.txt")
        table3 = np.array(table3)

        min_radius = self.earth_radius-max(table2[:,0])/1e3
        max_radius = self.earth_radius-min(table2[:,0])/1e3
        
        table=np.array(filter(lambda x: (x[0]>=min_radius and x[0]<=max_radius), table))

        assert(len(table) == len(table2))
        assert(len(table) == len(table3))

        self.table_radius = table[:,0]
        self.table_pressure = table[:,1] * 1e9
        self.table_density = table[:,2]
        self.table_vp = table3[:,1] / 1e3
        self.table_vs = table2[:,1] / 1e3



# this uses prem_withQ table
class prem_test(radiustable):
    def __init__(self):
        radiustable.__init__(self)

        table = tools.read_table("data/prem_withQ.txt")#data is: radius pressure density V_p V_s Q_K Q_mu
        table = np.array(table)
        self.table_radius = table[:,0] / 1e3 #convert to km 
        self.table_pressure = table[:,1] * 1e9
        self.table_density = table[:,2]
        self.table_vp = table[:,3] / 1e3
        self.table_vs = table[:,4] / 1e3


    


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
    p = np.arange(1e9,360e9,5e9)
    depths = map(s.depth, p) 
    #we could also just specify some depth levels directly like this:
    #depths = np.arange(35,5600,100)

    #now evaluate everything at the given depths levels (using interpolation)
    pressures, density, v_p, v_s, v_phi = s.evaluate_all_at(depths)

    # plot vs and vp and v_phi (note that v_phi is computed!)
    plt.plot(depths,v_p,'+-r', label='v_p')
    plt.plot(depths,v_s,'+-b', label='v_s')
    plt.plot(depths,v_phi,'--g', label='v_phi')
    plt.legend()
    plt.xlabel('depth in km')
    plt.ylabel('km/s')
    plt.show()

    s1=prem()
    depths=s1.internal_depth_list()
    pressures, density, v_p, v_s, v_phi = s1.evaluate_all_at(depths)
    plt.plot(depths,v_p,'+-r', label='v_p')
    plt.plot(depths,v_s,'+-b', label='v_s')
    plt.plot(depths,v_phi,'--g', label='v_phi')

    s2=prem_test()
    depths=s2.internal_depth_list()
    pressures, density, v_p, v_s, v_phi = s2.evaluate_all_at(depths)
    plt.plot(depths,v_p,'x-r', label='v_p')
    plt.plot(depths,v_s,'x-b', label='v_s')
    plt.plot(depths,v_phi,'x-g', label='v_phi')
    plt.show()
    
    print min(depths),max(depths)
