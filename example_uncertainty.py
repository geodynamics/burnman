# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""

"""

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..')) 

import burnman
from burnman import minerals
import misc.colors as colors 


if __name__ == "__main__":    
    figure=plt.figure(dpi=100,figsize=(12,10))
    prop={'size':12}
    dashstyle2=(6,3)
    dashstyle3=(10,2,2,2)
    dashstyle4=(4,9)

    seismic_model = burnman.seismic.prem() # pick from .prem() .slow() .fast() (see burnman/seismic.py)
    number_of_points = 10 #set on how many depth slices the computations should be done
    depths = np.linspace(850e3,2700e3, number_of_points)
    seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)

    

    def eval(uncertain):
        rock = burnman.composite([(minerals.SLB_2011.stishovite_err(uncertain), 1.0)] )
        rock.set_method('slb3')

        temperature = burnman.geotherm.adiabatic(seis_p,1900,rock)
        
        mat_rho, mat_vp, mat_vs, mat_vphi, mat_K, mat_G = \
            burnman.velocities_from_rock(rock, seis_p, temperature, burnman.averaging_schemes.voigt_reuss_hill())

        return seis_p, mat_vs, mat_vphi, mat_rho

    len = 9
    
    p, base_vs, base_vphi, _ = eval(np.zeros(len))



    names = ['$K_0$', '$K_0\'$', '$G_0$', '$G_0\'$','$\\theta_0$','$\gamma_0$','$q_0$','$\eta_0s$']

  
    for i in range(0,8):

        vsmin = base_vs
        vsmax = base_vs
        vphimin = base_vphi
        vphimax = base_vphi

        testrange = [-1.0, -0.5, 0.5, 1.0] #this seems to be enough for now
        for x in testrange:
            print i, names[i], x
            uncertain = np.zeros(len)
            uncertain[i]=x
            _, vs, vphi, _ = eval(uncertain)
            vsmin = np.minimum(vs,vsmin)
            vsmax = np.maximum(vs,vsmax)
            vphimin = np.minimum(vphi,vphimin)
            vphimax = np.maximum(vphi,vphimax)

        ax = figure.add_subplot(3,3,i+1)
        #plt.subplots_adjust(left=0, bottom=None, right=0, top=None, wspace=None, hspace=None)
        plt.subplots_adjust(wspace=0., hspace=0.2)
        
	print np.min((vsmin-base_vs)/base_vs)
        print np.max((vsmax-base_vs)/base_vs)
        print np.min((vphimin-base_vphi)/base_vphi)
        print np.max((vphimax-base_vphi)/base_vphi)


#        plt.plot(seis_p/1.e9,seis_vs/1.e3,color='k',dashes=dashstyle4,linewidth=1.0,marker='o', markersize=6,markerfacecolor='None',label='PREM')

        plt.plot(seis_p/1.e9,base_vs/1.e3,color=colors.color(3),linestyle=":",linewidth=1.0, markersize=6,markerfacecolor='None',label='pv')

        plt.plot(seis_p/1.e9,vsmin/1.e3,color=colors.color(3),dashes=dashstyle2,linewidth=1.0, markersize=6,markerfacecolor='None',label='min')
        plt.plot(seis_p/1.e9,vsmax/1.e3,color=colors.color(3),dashes=dashstyle2,linewidth=1.0, markersize=6,markerfacecolor='None',label='max')
#        plt.title('Vs %s +/- %d\\%% '%(names[i], spread[i]*100) )
        #plt.ylim([6.2,7.6])

        if (i%3==0):
            plt.ylabel('Wave speed (km/s)')
        else:
            ax.yaxis.set_ticklabels([])

        if (i>5):
            plt.xlabel('Pressure (GPa)')
        else:
            ax.xaxis.set_ticklabels([])

        #plt.subplot(3,3,i+1)#+10
#        plt.plot(seis_p/1.e9,seis_vphi/1.e3,color='k',dashes=dashstyle4,linewidth=1.0,marker='o', markersize=6,markerfacecolor='None',label='PREM')
        plt.plot(seis_p/1.e9,base_vphi/1.e3,color=colors.color(1),linestyle=':',linewidth=1.0, markersize=6,markerfacecolor='None',label='pv')
        plt.plot(seis_p/1.e9,vphimin/1.e3,color=colors.color(1),linestyle='-',linewidth=1.0, markersize=6,markerfacecolor='None',label='min')
        plt.plot(seis_p/1.e9,vphimax/1.e3,color=colors.color(1),linestyle='-',linewidth=1.0, markersize=6,markerfacecolor='None',label='max')
    
        plt.title(names[i] )
        #plt.ylim([8.5,12.])

        #plt.ylim([6.1,11.8])
        #plt.xlim([30,130])

#        if (reorder[i]==8):
#            handles, labels = ax.get_legend_handles_labels()
#            plt.legend((handles[0],handles[6],handles[2]), ['PREM','$V_\phi$','$V_S$'], loc='center right',prop=prop)



#    plt.savefig("uncertain.pdf",bbox_inches='tight')
    plt.show()


