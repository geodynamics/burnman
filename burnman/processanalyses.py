# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2019 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import

import warnings
import numpy as np
from scipy.optimize import curve_fit, minimize
from collections import Counter

from . import SolidSolution

def fit_composition(fitted_elements, composition, compositional_uncertainties, formulae, endmember_site_occupancies, normalize=True, name=None):
    """
    It is assumed that any elements not in composition were not measured (but may exist in unknown quantities).
    If distinct oxidation states or site occupancies were measured
    (by Moessbauer, for example), then formulae should be modified 

    fitted elements should be a list of strings
    composition and compositional uncertainties should be given as arrays
    formulae should either be given as arrays (with columns corresponding to elements) or as dictionaries
    If the compositional uncertainties can either be sigmas or covariance matrices

    The composition and associated uncertainties in endmember *amounts*, not proportions. If normalize=True, then the endmember amounts are normalized to a total of one.
    """

    if type(formulae[0]) is dict or type(formulae[0]) is Counter:
        stoichiometric_matrix = np.array([[f[e] if e in f else 0. for e in fitted_elements] for f in formulae])
    else:
        stoichiometric_matrix = formulae

    b = composition
    
    if len(compositional_uncertainties.shape) == 1:
        b_uncertainties = np.diag(compositional_uncertainties *
                                  compositional_uncertainties)
    else:
        b_uncertainties = compositional_uncertainties

    if np.linalg.det(b_uncertainties) < 1.e-30: # ensure uncertainty matrix is not singular
        #warnings.warn('The compositional covariance matrix for the {0} solution is nearly singular or not positive-definite (determinant = {1}). '
        #              'This is likely to be because your fitting parameters are not independent. '
        #              'For now, we increase all diagonal components by 1%. '
        #              'However, you may wish to redefine your problem.'.format(name, np.linalg.det(b_uncertainties)))

        b_uncertainties += np.diag(np.diag(b_uncertainties))*0.01

    endmember_constraints = lambda site_occ: [{'type': 'ineq', 'fun': lambda x, eq=eq: eq.dot(x)}
                                              for eq in site_occ]
    cons = endmember_constraints(endmember_site_occupancies.T)    
    A = stoichiometric_matrix.T
    
    fn = lambda A, *proportions: A.dot(proportions)
    popt, pcov = curve_fit(fn, A, b,
                           p0=np.array([0. for i in range(len(A.T))]),
                           sigma=b_uncertainties, absolute_sigma=True)

    res = np.sqrt((A.dot(popt) - b).dot(np.linalg.solve(b_uncertainties, A.dot(popt) - b)))

    # Check constraints
    if any([c['fun'](popt)<0. for c in cons]):
        warnings.warn('Warning: Simple least squares predicts an unfeasible solution composition for {0} solution. '
                      'Recalculating with site constraints. The covariance matrix must be treated with caution.'.format(name))
        fn = lambda x, A, b, b_uncertainties: np.sqrt((A.dot(popt) - b).dot(np.linalg.solve(b_uncertainties, A.dot(popt) - b)))
        sol = minimize(fn, popt, args=(A, b, b_uncertainties), method='SLSQP',constraints=cons)
        popt = sol.x
        res = sol.fun

    if normalize:
        sump = sum(popt)
        popt /= sump
        pcov /= sump*sump
        res /= sump
    return (popt, pcov, res)



def compute_and_set_phase_compositions(assemblage, verbose=False):
    for phase in assemblage.phases:
        if isinstance(phase, SolidSolution):
            try:
                name = phase.name
            except:
                name = 'unnamed solution' 
            popt, pcov, res = fit_composition(phase.fitted_elements,
                                              phase.composition,
                                              phase.compositional_uncertainties,
                                              phase.endmember_formulae,
                                              phase.solution_model.endmember_occupancies,
                                              normalize=True, name=name)
            
            # Convert uncertainties in amounts into uncertainties in proportions
            n_mbrs = len(phase.endmember_formulae)
            dpdx = np.zeros((n_mbrs, n_mbrs))
            for i in range(n_mbrs):
                for j in range(n_mbrs):
                    dpdx[i,j] = 1. - popt[i]/sum(popt) if i==j else -popt[i]/sum(popt)
            Cov_p = dpdx.dot(pcov).dot(dpdx.T)

            phase.set_composition(popt)
            phase.molar_fraction_covariances = Cov_p

            if verbose:
                print(phase.name)
                for i in range(n_mbrs):
                    print('{0}: {1:.3f} +/- {2:.3f}'.format(phase.endmember_names[i], popt[i], np.sqrt(Cov_p[i][i])))

def assemblage_affinity_misfit(assemblage):
    reaction_matrix = assemblage.stoichiometric_matrix(calculate_subspaces=True)[1]
    if len(reaction_matrix) == 0:
        print('No reactions between the phases'
              ' in this assemblage: {0}'.format([ph.name for ph in assemblage.phases]))
        return 0.
    else:
        # d(partial_gibbs)i/d(variables)j can be split into blocks
        n_mbrs = sum([phase.n_endmembers if isinstance(phase, SolidSolution) else 1 for phase in assemblage.phases])
        Cov_mu = np.zeros((n_mbrs, n_mbrs))
        dmudPT = np.zeros((n_mbrs, 2))
        mu = np.zeros(n_mbrs)
        
        i=0
        for phase in assemblage.phases:
            if isinstance(phase, SolidSolution):
                Cov_mu[i:i+phase.n_endmembers,i:i+phase.n_endmembers] = (phase.gibbs_hessian).dot(phase.molar_fraction_covariances).dot(phase.gibbs_hessian.T) # the hessian is symmetric, so transpose only taken for legibility...
                dmudPT[i:i+phase.n_endmembers] = np.array([phase.partial_volumes, -phase.partial_entropies]).T
                mu[i:i+phase.n_endmembers] = phase.partial_gibbs
                i += phase.n_endmembers
            else:
                dmudPT[i] = np.array([phase.V, -phase.S])
                mu[i] = phase.gibbs
                i += 1

        Cov_mu += dmudPT.dot(assemblage.state_covariances).dot(dmudPT.T)

        # Finally, we use the reaction matrix
        # (the nullspace of the stoichiometric matrix)
        # to calculate the affinities
        a = reaction_matrix.dot(mu)
        Cov_a = reaction_matrix.dot(Cov_mu).dot(reaction_matrix.T)
        try:
            chi_sqr = a.dot(np.linalg.solve(Cov_a, a))
        except:
            print('Could not find misfit for assemblage')
            try:
                print(assemblage.experiment_id)
            except:
                pass
            print([ph.name for ph in assemblage.phases])
            print(Cov_a)
            exit()
    return chi_sqr
