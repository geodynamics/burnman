# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import tools
import scipy.interpolate
import scipy.integrate
import matplotlib.pyplot as plt

class Seismic1DModel(object):
    """
    Base class for all the seismological models.
    """
    def __init__(self):
        pass



    def evaluate_all_at(self, depth_list):
        """
        Returns the lists of data for a Seismic1DModel for the depths provided
            
        Parameters
        ----------
        depth_list : array of floats
            Array of depths [m] to evaluate seismic model at.
            
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
    


    def internal_depth_list(self,mindepth=0.,maxdepth=6371.e3):
        """
        Returns a sorted list of depths where this seismic data is specified at. This allows you to compare the seismic data without interpolation.
        
        
        Returns
        -------
        depths : array of floats
        Depths [m].
        """
        raise ValueError, "not implemented"
        return 0

    def pressure(self, depth):
        """
        Parameters
        ----------
        depth : float or array of floats
        Depth(s) [m] to evaluate seismic model at.
        
        Returns
        -------
        pressure : float or array of floats
        Pressure(s) at given depth(s) in [Pa].
        """
        raise ValueError, "not implemented"
        return 0

    def v_p(self, depth):
        """
        Parameters
        ----------
        depth : float or array of floats
            Depth(s) [m] to evaluate seismic model at.
        
        Returns
        -------
        v_p : float or array of floats
            P wave velocity at given depth(s) in [m/s].
        """
        raise ValueError, "not implemented"
        return 0

    def v_s(self, depth):
        """
        Parameters
        ----------
        depth : float or array of floats
        Depth(s) [m] to evaluate seismic model at.
        
        Returns
        -------
        v_s : float or array of floats
        S wave velocity at given depth(s) in [m/s].
        """
        raise ValueError, "not implemented"
        return 0
    
    def v_phi(self, depth):
        """
        Parameters
        ----------
        depth_list : float or array of floats
        Depth(s) [m] to evaluate seismic model at.
        
        Returns
        -------
        v_phi : float or array of floats
        bulk sound wave velocity at given depth(s) in [m/s].
        """
        v_s=self.v_s(depth)
        v_p=self.v_p(depth)
        return np.sqrt(v_p*v_p-4./3.*v_s*v_s)
        
        
    def density(self, depth):
        """
        Parameters
        ----------
        depth : float or array of floats
        Depth(s) [m] to evaluate seismic model at.
        
        Returns
        -------
        density : float or array of floats
        Density at given depth(s) in [kg/m^3].
        """
        raise ValueError, "not implemented"
        return 0
    
    
    def G(self, depth):
        """
        Parameters
        ----------
        depth : float or array of floats
        Shear modulus at given for depth(s) in [Pa].
        """
        return np.power(self.v_s(depth),2.) * self.density(depth)
    
    def K(self, depth):
        """
        Parameters
        ----------
        depth : float or array of floats
        Bulk modulus at given for depth(s) in [Pa]
        """
        return np.power(self.v_phi(depth),2.) * self.density(depth)


    def QK(self, depth):
        """
        Parameters
        ----------
        depth : float or array of floats
        Depth(s) [m] to evaluate seismic model at.

        Returns
        -------
        Qk : float or array of floats
        Quality factor (dimensionless) for bulk modulus at given depth(s).
        """
        raise ValueError, "not implemented"
        return 0

    def QG(self, depth):
        """
        Parameters
        ----------
        depth : float or array of floats
        Depth(s) [m] to evaluate seismic model at.
        
        Returns
        -------
        QG : float or array of floats
        Quality factor (dimensionless) for shear modulus at given depth(s).
        """
        raise ValueError, "not implemented"
        return 0



    def depth(self, pressure):
        """
        Parameters
        ----------
        pressure : float or array of floats
        Pressure(s) [Pa] to evaluate depth at.

        Returns
        -------
        depth : float or array of floats
        Depth(s) [m] for given pressure(s)
        """
        raise ValueError, "not implemented"
        return -1
    


    def gravity (self, depth):
        """
        Parameters
        ----------
        depth : float or array of floats
        Depth(s) [m] to evaluate gravity at.
        
        Returns
        -------
        gravity : float or array of floats
        Gravity for given depths in [m/s^2]
        """
        raise ValueError, "not implemented"
        return -1


class SeismicTable(Seismic1DModel):
    """ 
    This is a base class that gets a 1D seismic model from a table indexed and
    sorted by radius. Fill the tables in the constructor after deriving
    from this class. This class uses :class:`burnman.seismic.Seismic1DModel`
    
    Note: all tables need to be sorted by increasing depth. self.table_depth needs to be defined
    Alternatively, you can also overwrite the _lookup function if you
    want to access with something else.
    """ 
    def __init__(self):
        Seismic1DModel.__init__(self)
        self.table_depth = []
        self.table_radius = []
        self.table_pressure = []
        #self.table_gravity = []
        self.table_density = []
        self.table_vp = []
        self.table_vs = []
        self.table_QG = []
        self.table_QK = []

        self.earth_radius = 6371.0e3
        
    def internal_depth_list(self,mindepth=0., maxdepth=6371.e3):
        return np.array([self.table_depth[x] for x in range(len(self.table_depth)) if self.table_depth[x]>mindepth and self.table_depth[x]<maxdepth])


    def pressure(self, depth):
        try:
            return self._lookup(depth, self.table_pressure)
        except:
            self._compute_pressure()
            return self._lookup(depth, self.table_pressure)

    def gravity(self, depth):
        try:
            return self._lookup(depth,self.table_gravity)
        except:
            self._compute_gravity()
            return self._lookup(depth,self.table_gravity)

    def v_p(self, depth):

        return self._lookup(depth, self.table_vp)

    def v_s(self, depth):

        return self._lookup(depth, self.table_vs)

    def QK(self, depth):

        return self._lookup(depth, self.table_QK)

    def QG(self, depth):

        return self._lookup(depth, self.table_QG)


    def density(self, depth):

        return self._lookup(depth, self.table_density)

    def depth(self, pressure):
        if pressure > max(self.table_pressure) or pressure < min(self.table_pressure)  :
            raise ValueError, "Pressure outside range of SeismicTable"

        depth = np.interp(pressure, self.table_pressure, self.table_depth )
        return depth

    def radius(self, pressure):

        radius = np.interp(pressure, self.table_pressure[::-1], self.earth_radius - self.table_depth[::-1])
        return radius
    
    def shift_discontinuities(self):
        #Shifting depths of discontinuities by 1 m, so they have unique depth values and interpolations work well....
        discontinuities=np.where(self.table_depth[1:]-self.table_depth[:-1]==0)[0]
        self.table_depth[discontinuities+1]=self.table_depth[discontinuities+1]+1.
        if not len(self.table_radius) == 0:
            self.table_radius[discontinuities+1]=self.table_radius[discontinuities+1]-1.


    def _lookup(self, depth, value_table):
        return np.interp(depth, self.table_depth, value_table)



    def _compute_gravity(self):
        #Calculate the gravity of the planet, based on a density profile.  This integrates
        #Poisson's equation in radius, under the assumption that the planet is laterally
        #homogeneous.
        #Create a spline fit of density as a function of radius
        density=self.table_density[::-1]
        radii=self.table_radius[::-1]
        rhofunc = scipy.interpolate.UnivariateSpline(radii, density)
 
        G = 6.67e-11
        #Numerically integrate Poisson's equation  --- can't get this part to work!!!!
        poisson = lambda p,x : 4.0 * np.pi * G * rhofunc(x) * x * x
        grav = np.ravel(scipy.integrate.odeint( poisson, 0.0, radii, full_output=0, printmessg=1))
        grav[1:] = grav[1:]/radii[1:]/radii[1:]
        grav[0] = 0.0 #Set it to zero a the center, since radius = 0 there we cannot divide by r^2
        
        
        #simple integration
        g=[]
        for i in range(len(radii)):
            g_tmp=scipy.integrate.trapz(G*4.*np.pi*density[0:i]*radii[0:i]*radii[0:i]/(radii[i])/(radii[i]),radii[0:i])
            g.append(g_tmp)
        self.table_gravity=g[::-1]


    def _compute_pressure(self):
        #Calculate the pressure profile based on density and gravity.  This integrates
        #the equation for hydrostatic equilibrium  P = rho g z.
        radii=self.table_radius
        density=self.table_density
        gravity=self.gravity(self.earth_radius-radii)
        #convert radii to depths
        depth = self.earth_radius-radii
        
        
        #### This isn't working and I'm not sure why...
        #Make a spline fit of density as a function of depth
        rhofunc = scipy.interpolate.UnivariateSpline( depth, density )
        #Make a spline fit of gravity as a function of depth
        gfunc = scipy.interpolate.UnivariateSpline( depth, gravity )
        #integrate the hydrostatic equation
        pressure = np.ravel(scipy.integrate.odeint( (lambda p, x : gfunc(x)* rhofunc(x)), 0.0,depth))
        
        
        p=[]
        for i in range(len(depth)):
            p_tmp=scipy.integrate.trapz(gravity[0:i]*density[0:i],depth[0:i])
            p.append(p_tmp)

        self.table_pressure=p


class PREM(SeismicTable):
    """
    Reads  PREM (1s) (input_seismic/prem.txt, :cite:`dziewonski1981`).
    See also :class:`burnman.seismic.SeismicTable`.
    """
    def __init__(self):
        SeismicTable.__init__(self)
        table = tools.read_table("input_seismic/prem.txt") # radius, pressure, density, v_p, v_s
        table = np.array(table)
        self.table_depth = table[:,0]
        self.table_radius = table[:,1]
        self.table_pressure = table[:,2]
        self.table_gravity = table[:,3]
        self.table_density = table[:,4]
        self.table_vp = table[:,5]
        self.table_vs = table[:,6]
        self.table_QK = table[:,7]
        self.table_QG = table[:,8]

        self.shift_discontinuities()



class Slow(SeismicTable):
    """
    Inserts the mean profiles for slower regions in the lower mantle (Lekic et al. 2012).
    We stitch together tables 'input_seismic/prem_lowermantle.txt', 'input_seismic/swave_slow.txt', 'input_seismic/pwave_slow.txt').
    See also :class:`burnman.seismic.SeismicTable`.
    """
    def __init__(self):
        SeismicTable.__init__(self)

        table = tools.read_table("input_seismic/prem.txt")#data is: radius pressure density V_p V_s Q_K Q_G
        table = np.array(table)
        table2 = tools.read_table("input_seismic/swave_slow.txt")
        table2 = np.array(table2)
        table3 = tools.read_table("input_seismic/pwave_slow.txt")
        table3 = np.array(table3)

        min_radius = self.earth_radius-max(table2[:,0])
        max_radius = self.earth_radius-min(table2[:,0])

        table=np.array(filter(lambda x: (x[1]>=min_radius and x[1]<=max_radius), table))


        self.table_depth = table[:,0]
        self.table_radius = table[:,1]
        self.table_pressure = table[:,2]
        self.table_density = table[:,4]
        self.table_vp = np.interp(self.table_depth,table3[:,0][::-1],table3[:,1][::-1])
        self.table_vs = np.interp(self.table_depth,table2[:,0][::-1],table2[:,1][::-1])

class Fast(SeismicTable):
    """
    Inserts the mean profiles for faster regions in the lower mantle (Lekic et al. 2012).
    We stitch together tables 'input_seismic/prem_lowermantle.txt', 'input_seismic/swave_fast.txt', 'input_seismic/pwave_fast.txt').
    See also :class:`burnman.seismic.Seismic1DModel`.
    """
    def __init__(self):
        SeismicTable.__init__(self)

        table = tools.read_table("input_seismic/prem.txt")#data is: radius pressure density V_p V_s Q_K Q_G
        table = np.array(table)
        table2 = tools.read_table("input_seismic/swave_fast.txt")
        table2 = np.array(table2)
        table3 = tools.read_table("input_seismic/pwave_fast.txt")
        table3 = np.array(table3)

        min_radius = self.earth_radius-max(table2[:,0])
        max_radius = self.earth_radius-min(table2[:,0])

        table=np.array(filter(lambda x: (x[1]>=min_radius and x[1]<=max_radius), table))


        self.table_depth = table[:,0]
        self.table_radius = table[:,1]
        self.table_pressure = table[:,2]
        self.table_density = table[:,4]
        self.table_vp = np.interp(self.table_depth,table3[:,0][::-1],table3[:,1][::-1])
        self.table_vs = np.interp(self.table_depth,table2[:,0][::-1],table2[:,1][::-1])


class REF(SeismicTable):
    """
        Reads  REF or STW05 (1s) (input_seismic/STW105.txt, :cite:`kustowski2008`).
        See also :class:`burnman.seismic.SeismicTable`.
        """
    def __init__(self):
        SeismicTable.__init__(self)
        table = tools.read_table("input_seismic/STW105.txt") # radius, pressure, density, v_p, v_s
        table = np.array(table)
        self.table_radius = table[:,0][::-1]
        self.table_density = table[:,1][::-1]
        self.table_vpv = table[:,2][::-1]
        self.table_vsv = table[:,3][::-1]
        self.table_QK = table[:,4][::-1]
        self.table_QG = table[:,5][::-1]
        self.table_vph = table[:,6][::-1]
        self.table_vsh = table[:,7][::-1]
        
        self.table_depth=self.earth_radius-self.table_radius
        
        # Voigt averages for Vs and Vp
        self.table_vs=np.sqrt((2.*self.table_vsv*self.table_vsv+self.table_vsh*self.table_vsh)/3.)
        self.table_vp=np.sqrt((self.table_vpv*self.table_vpv+4.*self.table_vph*self.table_vph)/5.)

        self.shift_discontinuities()


class IASP91(SeismicTable):
    """
        Reads  REF/STW05 (input_seismic/STW105.txt, :cite:`kustowski2008`).
        See also :class:`burnman.seismic.SeismicTable`.
        """
    def __init__(self):
        SeismicTable.__init__(self)
        table = tools.read_table("input_seismic/IASP91.txt") # radius, pressure, density, v_p, v_s
        table = np.array(table)
        self.table_depth = table[:,0]
        self.table_radius = table[:,1]
        self.table_vp = table[:,2]
        self.table_vs = table[:,3]
        
        
        self.shift_discontinuities()

class AK135(SeismicTable):
    """
        Reads  AK135 (input_seismic/ak135.txt, :cite:`kennett1995`).
        See also :class:`burnman.seismic.SeismicTable`.
        """
    def __init__(self):
        SeismicTable.__init__(self)
        table = tools.read_table("input_seismic/ak135.txt") # radius, pressure, density, v_p, v_s
        table = np.array(table)
        self.table_depth = table[:,0]
        self.table_radius = table[:,1]
        self.table_density=table[:,2]
        self.table_vp = table[:,3]
        self.table_vs = table[:,4]
        self.table_QG = table[:,5]
        self.table_QK = table[:,6]
        
        self.shift_discontinuities()


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
    v_p *= 1. - 1. / 2. * cot * 1. / Qp  # Matas et al. (2007) page 1
    v_s *= 1. - 1. / 2. * cot * 1. / Qs
    v_phi *= 1. - 1. / 2. * cot * 1. / Qphi
    return v_p, v_s, v_phi

"""
shared variable of prem, so that other routines do not need to create
prem over and over. See geotherm for example.
"""
prem_model = PREM()

