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
import colors

class my_perovskite(burnman.material):
    """
    based on Stixrude & Lithgow-Bertelloni 2011 and references therein  
    """
    def __init__(self, uncertain):
        self.params = {
            'equation_of_state':'slb3',
            'ref_V': 24.45e-6,
            'ref_K': 251.0e9 * uncertain[0],   
            'K_prime': 4.1 * uncertain[1],
            'ref_G': 173.0e9 * uncertain[2],
            'G_prime': 1.7 * uncertain[3],
            'molar_mass': .1000,
            'n': 5,
            'ref_Debye': 905. * uncertain[4], # less important?
            'ref_grueneisen': 1.57 * uncertain[5],
            'q0': 1.1 * uncertain[6],
            'eta_0s': 2.6 * uncertain[7]}

if __name__ == "__main__":    
    figure=plt.figure(dpi=100,figsize=(12,10))
    prop={'size':12}
    plt.rc('text', usetex=True)
    plt.rcParams['text.latex.preamble'] = '\usepackage{relsize}'
    plt.rc('font', family='sanserif')

    dashstyle2=(6,3)
    dashstyle3=(10,2,2,2)
    dashstyle4=(4,9)

    seismic_model = burnman.seismic.prem() # pick from .prem() .slow() .fast() (see burnman/seismic.py)
    number_of_points = 10 #set on how many depth slices the computations should be done
    depths = np.linspace(850e3,2700e3, number_of_points)
    seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)

    

    def eval(uncertain):
        rock = burnman.composite ( [ (my_perovskite(uncertain), 1.0) ])
        rock.set_method('slb3')

        temperature = burnman.geotherm.adiabatic(seis_p,1900*uncertain[8],rock)
        
        mat_rho, mat_vp, mat_vs, mat_vphi, mat_K, mat_G = \
            burnman.velocities_from_rock(rock, seis_p, temperature, burnman.averaging_schemes.voigt_reuss_hill())

        return seis_p, mat_vs, mat_vphi, mat_rho

    len = 9
    
    p, base_vs, base_vphi, _ = eval(np.ones(len))

    spread = [.1, .1, .1, .1, .1, .1, .1, .1, .1]


    names = ['$K_0$', '$K\'$', '$G_0$', '$G\'$', '$\\theta_0$', '$\gamma_0$', '$q_0$', '$\eta_{s0}$','$T_0$']

    reorder = [0,1,3,4,5,6,7,8,2]

    for i in range(0,len):

        vsmin = base_vs
        vsmax = base_vs
        vphimin = base_vphi
        vphimax = base_vphi

        #testrange = [-1, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 1.0]
        testrange = [-1.0, -0.5, -0.1, 0.1, 0.5, 1.0] #this seems to be enough for now
        for x in testrange:
            print i, names[i], x
            uncertain = np.ones(len)
            uncertain[i]=uncertain[i]+spread[i]*x
            _, vs, vphi, _ = eval(uncertain)
            vsmin = np.minimum(vs,vsmin)
            vsmax = np.maximum(vs,vsmax)
            vphimin = np.minimum(vphi,vphimin)
            vphimax = np.maximum(vphi,vphimax)

        ax = figure.add_subplot(3,3,reorder[i]+1)
        #plt.subplots_adjust(left=0, bottom=None, right=0, top=None, wspace=None, hspace=None)
        plt.subplots_adjust(wspace=0, hspace=0.2)
        

        plt.plot(seis_p/1.e9,seis_vs/1.e3,color='k',dashes=dashstyle4,linewidth=1.0,marker='o', markersize=6,markerfacecolor='None',label='PREM')

        plt.plot(seis_p/1.e9,base_vs/1.e3,color=colors.color(3),linestyle=":",linewidth=1.0, markersize=6,markerfacecolor='None',label='pv')

        plt.plot(seis_p/1.e9,vsmin/1.e3,color=colors.color(3),dashes=dashstyle2,linewidth=1.0, markersize=6,markerfacecolor='None',label='min')
        plt.plot(seis_p/1.e9,vsmax/1.e3,color=colors.color(3),dashes=dashstyle2,linewidth=1.0, markersize=6,markerfacecolor='None',label='max')
        plt.title('Vs %s +/- %d\\%% '%(names[i], spread[i]*100) )
        plt.ylim([6.2,7.6])

        if (reorder[i]%3==0):
            plt.ylabel('km/s')
        else:
            ax.yaxis.set_ticklabels([])

        if (reorder[i]>5):
            plt.xlabel('Pressure (GPa)')
        else:
            ax.xaxis.set_ticklabels([])

        #plt.subplot(3,3,reorder[i]+1)#+10
        plt.plot(seis_p/1.e9,seis_vphi/1.e3,color='k',dashes=dashstyle4,linewidth=1.0,marker='o', markersize=6,markerfacecolor='None',label='PREM')
        plt.plot(seis_p/1.e9,base_vphi/1.e3,color=colors.color(1),linestyle=':',linewidth=1.0, markersize=6,markerfacecolor='None',label='pv')
        plt.plot(seis_p/1.e9,vphimin/1.e3,color=colors.color(1),linestyle='-',linewidth=1.0, markersize=6,markerfacecolor='None',label='min')
        plt.plot(seis_p/1.e9,vphimax/1.e3,color=colors.color(1),linestyle='-',linewidth=1.0, markersize=6,markerfacecolor='None',label='max')

        plt.title('%s $\pm %d\\%%$ '%(names[i], spread[i]*100) )
        #plt.ylim([8.5,12.])

        plt.ylim([6.1,11.8])
        plt.xlim([30,130])

        if (reorder[i]==8):
            handles, labels = ax.get_legend_handles_labels()
            plt.legend((handles[0],handles[6],handles[2]), ['PREM','$V_\phi$','$V_S$'], loc='center right',prop=prop)



    plt.savefig("uncertain.pdf",bbox_inches='tight')
    #plt.show()


