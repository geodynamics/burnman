# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
example_seismic
---------------

Shows the various ways to input seismic models (:math:`V_s, V_p, V_{\phi}, \\rho`) as a
function of depth (or pressure) as well as different velocity model libraries
available within Burnman:

1. PREM :cite:`dziewonski1981`
2. Reference model for fast regions (outside the LLSVP's) in the lower mantle
   :cite:`Lekic2012`
3. Reference model for slow regions (LLSVP's) in the lower mantle :cite:`Lekic2012`

This example will first calculate or read in a seismic model and plot the
model along the defined pressure range. The example also illustrates how to import a seismic model of your choice, here shown by importing AK135 :cite:`Kennett1995`.

*Uses:*

* :doc:`seismic`



*Demonstrates:*

* Utilization of library seismic models within BurnMan
* Input of user-defined seismic models


"""

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..'))

import burnman

if __name__ == "__main__":

    #create a seismic dataset from prem:
    s=burnman.seismic.PREM()

    # specify where we want to evaluate, here we map from pressure to depth
    #format p = np.arange (starting pressure, ending pressure, pressure step) (in Pa)
    p = np.arange(1.0e9,360.0e9,5e9)
    depths = np.array([s.depth(pr) for pr in p])
    #we could also just specify some depth levels directly like this:
    #depths = np.arange(35e3,5600e3,100e3)
    #we could also use the data points where the seismic model is specified:
    depths = s.internal_depth_list()

    #now evaluate everything at the given depths levels (using interpolation)
    pressures, density, v_p, v_s, v_phi = s.evaluate_all_at(depths)
    
    # plot vs and vp and v_phi (note that v_phi is computed!)
    plt.subplot(2,2,1)
    plt.title('prem')
    plt.plot(depths/1.e3,v_p/1.e3,'+-r', label='v_p')
    plt.plot(depths/1.e3,v_s/1.e3,'+-b', label='v_s')
    plt.plot(depths/1.e3,v_phi/1.e3,'--g', label='v_phi')
    plt.legend(loc='lower left')
    plt.xlabel('depth in km')
    plt.ylabel('km/s')

    # plot pressure,density vs depth from prem:
    plt.subplot(2,2,2)
    plt.title('prem')
    plt.plot(depths/1.e3,pressures/1.e9,'-r', label='pressure')
    plt.ylabel('GPa')
    plt.xlabel('depth in km')
    plt.legend(loc='upper left')
    plt.twinx()
    plt.ylabel('g/cc')
    plt.plot(depths/1.e3,density/1.e3,'-b', label='density')
    plt.legend(loc='lower right')


    #now load fast and slow regionalized models (Lekic et al. 2012):
    sslow = burnman.seismic.Slow()
    depths2 = sslow.internal_depth_list()
    pressures2, density2, v_p2, v_s2, v_phi2 = sslow.evaluate_all_at(depths2)

    sfast = burnman.seismic.Fast()
    depths3 = sfast.internal_depth_list()
    pressures3, density3, v_p3, v_s3, v_phi3 = sfast.evaluate_all_at(depths3)


    # plotting
    plt.subplot(2,2,3)
    plt.plot(pressures/1.e9,v_p/1.e3,'-k', label='v_p prem')
    plt.plot(pressures2/1.e9,v_p2/1.e3,'-r', label='v_p slow')
    plt.plot(pressures3/1.e9,v_p3/1.e3,'-b', label='v_p fast')

    plt.legend(loc='lower right')
    plt.xlim([30,136])
    plt.ylim([11,14])
    plt.xlabel('pressure')
    plt.ylabel('km/s')

    plt.subplot(2,2,4)
    plt.plot(pressures/1.e9,v_s/1.e3,'-k', label='v_s prem')
    plt.plot(pressures2/1.e9,v_s2/1.e3,'-r', label='v_s slow')
    plt.plot(pressures3/1.e9,v_s3/1.e3,'-b', label='v_s fast')

    plt.legend(loc='upper left')
    plt.xlim([30,136])
    plt.ylim([6,8])
    plt.xlabel('pressure')
    plt.ylabel('km/s')

    plt.show()

    # Loading an a seismic model from a file. In this case AK135 (Kennett et al. 1995).
    # Note that file input is assumed to be in SI units
    plt.close()
    # Load data table

    class ak135_table(burnman.seismic.SeismicRadiusTable):
        def __init__(self):
            burnman.seismic.SeismicRadiusTable.__init__(self)
            # In format: radius, pressure, density, v_p, v_s
            table = burnman.tools.read_table("input_seismic/ak135_lowermantle.txt")
            table = np.array(table)
            self.table_radius = table[:,0]
            self.table_pressure = table[:,1]
            self.table_density = table[:,2]
            self.table_vp = table[:,3]
            self.table_vs = table[:,4]


    ak=ak135_table()
    # specify where we want to evaluate, here we map from pressure to depth
    depths = np.linspace(700e3, 2800e3, 40)
    #now evaluate everything at the given depths levels (using interpolation)
    pressures, density, v_p, v_s, v_phi = ak.evaluate_all_at(depths)
    # plot vs and vp and v_phi (note that v_phi is computed!)
    plt.subplot(2,2,1)
    plt.title('ak135')
    plt.plot(depths/1.e3,v_p/1.e3,'+-r', label='v_p')
    plt.plot(depths/1.e3,v_s/1.e3,'+-b', label='v_s')
    plt.plot(depths/1.e3,v_phi/1.e3,'--g', label='v_phi')
    plt.legend(loc='lower left')
    plt.xlabel('depth in km')
    plt.ylabel('km/s')

    # plot pressure,density vs depth from prem:
    plt.subplot(2,2,2)
    plt.title('ak135')
    plt.plot(depths/1.e3,pressures/1.e9,'-r', label='pressure')
    plt.ylabel('GPa')
    plt.xlabel('depth in km')
    plt.legend(loc='upper left')
    plt.twinx()
    plt.ylabel('g/cc')
    plt.plot(depths/1.e3,density/1.e3,'-b', label='density')
    plt.legend(loc='lower right')
    plt.savefig("output_figures/example_seismic.png")
    plt.show()
