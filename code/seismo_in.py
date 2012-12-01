# INPUT for the seismic velocity and density model
# Density will always be PREM for now
#standard numpy, scipy imports
import numpy as np
import matplotlib.pyplot as plt

#eos imports
import seismic_models as seis

class seismo_in:
    def __init__(self, ):
        self.params = {	'name':'prem',  ## choose model: 'prem','ref_fast','ref_slow'
                        'attenuation_correction':'off', ## if 'on' this (small) correction will correct for the value in QL6 (Durek & Ekstrom) assuming a Poisson solid for long period (10^2s) waves (see page 4 in Matas et al. 2007). Attenuation corrections are arguably needed when comparing measurements done at GHz to seismic observations on the order of mHz-Hz. However effects of attuation are ignorable above 1Hz. If you choose 'prem', which is the version of PREM for 1Hz, then no correction is needed.
                        'depth_min':0.0, # minimum depth considered in km, choose 0.0 to use the limits given by the seismic model
                        'depth_max':0.0, # maximum depth considered in km, choose 0.0 to use the limits given by the seismic model
                        'depth_step':100.0, # steps at which the seismic velocity and calculated velocity are compared (in between is interpolated)
            }

    def depth_bounds(self):
	global depth_bounds_ran_before
        [dmin,dmax]=seis.lookup_depthbounds(self.params['name'])
        dmin=dmin/1e3;
        dmax=dmax/1e3;
	if (self.params['depth_max']==0.0):
            depth_max=dmax
            depth_min=dmin
        elif (dmin>self.params['depth_min'] or dmax<self.params['depth_max']):    
                raise ValueError("your chosen seismic model is bounded at",dmin," km and",dmax," km. Set depthbounds to 0.0 in seismo_in to use the bounds set by the model")
        else:
            	depth_max=self.params['depth_max']
            	depth_min=self.params['depth_min']
	try: depth_bounds_ran_before
        except:
            print "minimum depth considered is", depth_min
            print "maximum depth considered is", depth_max
            depth_bounds_ran_before=True
        return depth_min, depth_max  
    def radius(self):
        global radius_ran_before
        radius=(np.arange(self.depth_bounds()[0],self.depth_bounds()[1],self.params['depth_step']))*1e3;
        try: radius_ran_before
	except:
	    print "calculations are done for", len(radius)," depths"
            radius_ran_before=True
        return radius    
    def pressure(self):
        r=self.radius()
        press=[seis.lookup_pressure(self.params['name'],r[y]) for y in range(len(r))]
        return press
    # put in dependence on seismic models, give in self                
    def density(self):
        p=self.pressure()
        rho=[seis.lookup_(self.params['name'],p[i],2) for i in range(len(p))]
        return rho            
    def v_s(self):
        p=self.pressure()  
        return [seis.lookup_(self.params['name'],p[i],4) for i in range(len(p))]                
    def v_p(self):
        p=self.pressure()
        return [seis.lookup_(self.params['name'],p[i],3) for i in range(len(p))]             
    def v_phi(self):
        vs=self.v_s()
        vp=self.v_p()
        return [np.sqrt(pow(vp[i],2.0)-4/3*pow(vs[i],2.0)) for i in range(len(vp))]







