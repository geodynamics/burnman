
# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
    
example_convert_eos_order.py
----------------

This example converts a mineral that has been fit with a second order EoS (slb2) to a third order equation of state



"""
import os, sys, numpy as np, matplotlib.pyplot as plt
import scipy.optimize as opt

#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..'))

import burnman

if __name__ == "__main__":
    #in form name_of_mineral (burnman.mineral <- creates list with parameters)
    class original_material(burnman.Mineral):
        def __init__(self):
            burnman.Mineral.__init__(self)
            self.params = {
                'equation_of_state': 'slb2', # original method
                'V_0': 24.45e-6,
                'K_0': 253.0e9,
                'Kprime_0': 4.1,
                'G_0': 172.9e9,
                'Gprime_0': 1.56,
                'molar_mass': .1000,
                'n': 5,
                'Debye_0': 1100.,
                'grueneisen_0': 1.4,
                'q_0': 1.4,
                'eta_s_0': 2.6 }
    

    def calc_velocities(mineral,method,pressures,temperature):
        
        mineral.set_method(method)
        
        mat_rho, mat_vp, mat_vs, mat_vphi, mat_K, mat_G = burnman.velocities_from_rock(mineral,pressures, temperature)
        return mat_K,  mat_rho, mat_G

    def error_K(guess, test_mineral, method,pressures, temperature, ori_K,ori_rho,ori_G):
        
        test_mineral.params['K_0']=guess[0]
        test_mineral.params['Kprime_0']=guess[1]

        K,rho,G = calc_velocities(test_mineral, method,pressures,temperature)

        K_err=burnman.l2(pressures,K,ori_K)
        G_err=burnman.l2(pressures,G,ori_G)
        rho_err=burnman.l2(pressures,rho,ori_rho)
        return K_err+rho_err

    def error_G(guess, test_mineral, method,pressures, temperature, ori_K,ori_rho,ori_G):
        

        test_mineral.params['G_0']=guess[0]
        test_mineral.params['Gprime_0']=guess[1]
        K,rho,G = calc_velocities(test_mineral, method,pressures,temperature)
        
        K_err=burnman.l2(pressures,K,ori_K)
        G_err=burnman.l2(pressures,G,ori_G)
        rho_err=burnman.l2(pressures,rho,ori_rho)
        return G_err


    pressures = np.linspace(25.e9, 135.e9, 10)
    temperature = burnman.geotherm.brown_shankland(pressures)
    # calculate moduli to fit
    ori_K, ori_rho, ori_G=calc_velocities(original_material(),'slb2',pressures,temperature)
    print ori_G
    min=original_material()
    #make the mineral to fit
    guess = [min.params['K_0'],min.params['Kprime_0']]
    print guess

    #first, do the second-order fit
    func = lambda x : error_K( x, min, 'slb3',pressures, temperature, ori_K, ori_rho, ori_G)
    sol_K = opt.fmin(func, guess)
    
    
    min.params['K_0']=sol_K[0]
    min.params['Kprime_0']=sol_K[1]
    guess = [min.params['G_0'],min.params['Gprime_0']]
    func = lambda x : error_G( x, min, 'slb3',pressures, temperature, ori_K, ori_rho, ori_G)
    sol_G = opt.fmin(func, guess)
    print sol_K,sol_G


    # calculate moduli with new solution
    minnew=original_material()
    minnew.params['K_0']=sol_K[0]
    minnew.params['Kprime_0']=sol_K[1]
    minnew.params['G_0']=sol_G[0]
    minnew.params['Gprime_0']=sol_G[1]
    new_K, new_rho, new_G=calc_velocities(original_material(),'slb3',pressures,temperature)

    plt.plot(pressures/1.e9,ori_K/1.e9,color='r', linestyle='-', linewidth=2, label = "K 2nd order")
    plt.plot(pressures/1.e9,new_K/1.e9,color='r', linestyle='-.', linewidth=2, label = "K 3rd order")
    plt.plot(pressures/1.e9,ori_G/1.e9,color='b', linestyle='-', linewidth=2, label = "G 2nd order")
    plt.plot(pressures/1.e9,new_G/1.e9,color='b', linestyle='-.', linewidth=2, label = "G 3rd order")
    #plt.ylim([6.55, 8])
    plt.xlim([25., 135.])
    plt.ylabel("Moduli (GPa)")
    plt.xlabel("Pressure (GPa)")
    plt.legend(loc = "lower right",prop={'size':12},frameon=False)
    plt.savefig("output_figures/example_fit_data.png")
    plt.show()
