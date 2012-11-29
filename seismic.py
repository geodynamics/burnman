import math

# Reuss-Voigt-Hill average
# inp: bulk modulus, shear modulus, density
# return K_Si and gamma_i 
def voigt_reuss_hill(molar_abundance, molar_weight, modulus, density, T):
    # from matas 2007, Appendix C

    n_phase = len(modulus)

    assert n_phase == len(molar_abundance)
    assert n_phase == len(molar_weight)
    assert n_phase == len(density)

    it = range(n_phase)

    n_i = molar_abundance  # molar abundance for phase i 
    M_i = molar_weight  # molar weight for phase i
    total_molar = sum(n_i)
    x_i = [(n_i[i] / total_molar) for i in it] # molar fraction of phase i
    
    #V_i = n_i * M_i / rho:
    V_i = [(n_i[i]*M_i[i]/density[i]) for i in it]
    #V = sum n_i V_i:
    V = sum(i*j for i, j in zip(n_i,V_i))
    #nu_i = x_i * V_i / V:
    nu_i = [(x_i[i] * V_i[i] / V) for i in it]

    #X_i = K_Si = K_Ti (1+alpha_i gamma_i T)
    #alpha_i: thermal expansion
    #gamma_i: grueneisen
    #K_Ti: bulk modulus    
    #not needed: X_i = [ (modulus[i] * (1.+thermal_exp[i] * grueneisen[i] * T)) for i in it]
    X_i = modulus

    #X_V = sum nu_i X_i:
    X_V = sum(i*j for i,j in zip(nu_i,X_i))
    #X_R = 1 / sum(nu_i/X_i):
    X_R = 1. / sum(i/j for i,j in zip(nu_i,X_i))

    return (X_V + X_R) / 2

# compute V_p, V_s, and V_phi using Voigt-Reuss-Hill
# input: molar_abundance, molar_weight, bulk_modulus, shear_modulus, density (arrays)
# input: T
# returns V_p,V_s,V_phi
def get_velocities(molar_abundance, molar_weight, bulk_modulus, shear_modulus, density, T):

    it = range(len(molar_abundance))
    n_i = molar_abundance  # molar abundance for phase i 
    M_i = molar_weight  # molar weight for phase i
    #V_i = n_i * M_i / density:
    V_i = [(n_i[i]*M_i[i]/density[i]) for i in it]
    #V = sum n_i V_i:
    V = sum(i*j for i, j in zip(n_i,V_i))    
    # avg_density = 1./ V sum(n_i M_i):
    avg_density = 1./ V * sum((n_i[i]*M_i[i]) for i in it)

    K_s = voigt_reuss_hill(molar_abundance, molar_weight, bulk_modulus, density, T)
    mu = voigt_reuss_hill(molar_abundance, molar_weight, shear_modulus, density, T)
    V_p = math.sqrt((K_s + 4./3. * mu) / avg_density )
    V_s = math.sqrt(mu / avg_density)
    V_phi = math.sqrt(K_s / avg_density)

    return V_p, V_s, V_phi
