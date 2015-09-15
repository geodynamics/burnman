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

    #List of seismic 1D models 
    models=[burnman.seismic.PREM(), burnman.seismic.STW105(),burnman.seismic.AK135(),burnman.seismic.IASP91()]
    colors = ['r','b','m','k']
    # Variables to plot
    vars=['pressure','gravity','v_p','v_s','v_phi','density']#,
    units=['Pa','m/s^2','m/s','m/s','m/s','kg/m^3','Pa','m/s^2']
    plt.figure(figsize=(14,8))
    

    # Run through models and variables
    for a in range(len(vars)):
        for m in range(len(models)):

                # specify where we want to evaluate, here we map from pressure to depth
                #1. format p = np.arange (starting pressure, ending pressure, pressure step) (in Pa)
                #p = np.arange(1.0e9,360.0e9,1.e9)
                #depths = np.array([models[m].depth(pr) for pr in p])
                #2. we could also just specify some depth levels directly like this
                #depths = np.arange(700e3,2800e3,100e3)
                #3. we could also use the data points where the seismic model is specified over a depth range, this will bring out any discontinuities
                #this is the preferred way to plot seismic discontinuities correctly
                depths = models[m].internal_depth_list(mindepth=-1.e3, maxdepth=6372.1e3)
                #now evaluate everything at the given depths levels (using interpolation)
                try:
                    values=getattr(models[m],vars[a])(depths)
                    #Plot
                    plt.plot(depths/1.e3,values,color=colors[m],linestyle='-',label=models[m].__class__.__name__)
                except:
                    print vars[a] + ' is not defined for ' + models[m].__class__.__name__

        plt.title(vars[a])
        plt.legend(loc='center right')
        plt.xlabel('depth in km')
        plt.ylabel(units[a])
        plt.show()


    
    ## Alternatively one is able to evaluate all the variables for a model in a single line
    #  values = model.evaluate(vars_list,depth_list)
    pressure, gravity, v_p, v_s, v_phi, density = models[0].evaluate(['pressure','gravity','v_p','v_s','v_phi','density'],models[0].internal_depth_list(mindepth=-1.e3, maxdepth=6372.1e3))





    ## The following shows how to read in your own model from a file
    ## Model needs to be defined with increasing depth and decreasing radius. In this case the table is switched.
    class ak135_table(burnman.seismic.SeismicTable):
        def __init__(self):
            burnman.seismic.SeismicTable.__init__(self)
            # In format: radius, pressure, density, v_p, v_s
            table = burnman.tools.read_table("input_seismic/ak135_lowermantle.txt")
            table = np.array(table)
            self.table_radius = table[:,0][::-1]
            self.table_pressure = table[:,1][::-1]
            self.table_density = table[:,2][::-1]
            self.table_vp = table[:,3][::-1]
            self.table_vs = table[:,4][::-1]

            #self.table_depth needs to be defined and needs to be increasing
            self.table_depth=self.earth_radius-self.table_radius


    ak=ak135_table()
    # specify where we want to evaluate, here we map from pressure to depth
    depths = np.linspace(700e3, 2800e3, 40)
    #now evaluate everything at the given depths levels (using interpolation)
    pressures, density, v_p, v_s, v_phi = ak.evaluate(['pressure','density','v_p','v_s','v_phi'],depths)
    # plot vs and vp and v_phi (note that v_phi is computed!)
    plt.figure(figsize=(10,5))
    plt.subplot(1,2,1)
    plt.title('ak135')
    plt.plot(depths/1.e3,v_p/1.e3,'+-r', label='v_p')
    plt.plot(depths/1.e3,v_s/1.e3,'+-b', label='v_s')
    plt.plot(depths/1.e3,v_phi/1.e3,'--g', label='v_phi')
    plt.legend(loc='lower left')
    plt.xlabel('depth in km')
    plt.ylabel('km/s')

    # plot pressure,density vs depth from prem:
    plt.subplot(1,2,2)
    plt.title('ak135')
    plt.plot(depths/1.e3,pressures/1.e9,'-r', label='pressure')
    plt.ylabel('GPa')
    plt.xlabel('depth in km')
    plt.legend(loc='upper left')
    plt.twinx()
    plt.ylabel('g/cc')
    plt.plot(depths/1.e3,density/1.e3,'-b', label='density')
    plt.legend(loc='lower right')
    plt.show()
