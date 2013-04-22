import os, sys, numpy as np, matplotlib.pyplot as plt

import scipy.optimize as opt
import burnman

#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
        sys.path.insert(1,os.path.abspath('..')) 

def calc_shear_velocities(ref_mu, mu_prime, mineral, pressures): 

    mineral.params['ref_mu'] = ref_mu
    mineral.params['mu_prime'] = mu_prime
    mineral.set_method('bm')

    shear_velocities = np.empty_like(pressures)
    for i in range(len(pressures)):
        mineral.set_state(pressures[i], 0.0) # set state with dummy temperature
        shear_velocities[i] = mineral.v_s()

    return shear_velocities

def error(guess, test_mineral, pressures, obs_vs):
    vs = calc_shear_velocities(guess[0], guess[1], test_mineral, pressures)

    vs_l2 = [ (vs[i] - obs_vs[i])*(vs[i] - obs_vs[i]) for i in range(len(obs_vs)) ]
    l2_error = sum(vs_l2)

    return l2_error


mg_perovskite_data = np.loadtxt("input_minphys/Murakami_perovskite.txt")
obs_pressures = mg_perovskite_data[:,0]*1.e9
obs_vs = mg_perovskite_data[:,2]*1000.

#make the mineral to fit
guess = [200.e9, 2.0]
mg_perovskite_test = burnman.material()
mg_perovskite_test.params['ref_V'] = 24.607e-6
mg_perovskite_test.params['ref_K'] = 251.9e9
mg_perovskite_test.params['K_prime'] = 4.01
mg_perovskite_test.params['molar_mass'] = .1053

func = lambda x : error( x, mg_perovskite_test, obs_pressures, obs_vs)
sol = opt.fmin(func, guess)

print sol 

pressures = np.linspace(30.e9, 130.e9, 100)
model_vs = calc_shear_velocities(sol[0], sol[1], mg_perovskite_test, pressures)

plt.plot(pressures/1.e9,model_vs/1000.,color='r')
plt.scatter(obs_pressures/1.e9, obs_vs/1000.)
plt.ylim([4, 8])
plt.title("Vs (km/s)")
plt.show()
