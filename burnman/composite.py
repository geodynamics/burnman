# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import warnings

from .material import Material, material_property
from .mineral import Mineral
from . import averaging_schemes
from . import chemicalpotentials
from .solidsolution import SolidSolution
from sympy import Matrix, nsimplify
import scipy.optimize as opt
from scipy.linalg import pinv

def check_pairs(phases, fractions):
        if len(fractions) < 1:
            raise Exception('ERROR: we need at least one phase')

        if len(phases) != len(fractions):
            raise Exception(
                'ERROR: different array lengths for phases and fractions')

        total = sum(fractions)
        if abs(total - 1.0) > 1e-10:
            raise Exception(
                'ERROR: list of molar fractions does not add up to one')
        for p in phases:
            if not isinstance(p, Mineral):
                raise Exception(
                    'ERROR: object of type ''%s'' is not of type Mineral' % (type(p)))


# static composite of minerals/composites
class Composite(Material):

    """
    Base class for a composite material.
    The static phases can be minerals or materials,
    meaning composite can be nested arbitrarily.

    The fractions of the phases can be input
    as either 'molar' or 'mass' during instantiation,
    and modified (or initialised) after this point by
    using set_fractions.

    This class is available as ``burnman.Composite``.
    """

    def __init__(self, phases, fractions=None, fraction_type='molar', name='Unnamed composite'):
        """
        Create a composite using a list of phases and their fractions (adding to 1.0).

        Parameters
        ----------
        phases: list of :class:`burnman.Material`
            list of phases.
        fractions: list of floats
            molar or mass fraction for each phase.
        fraction_type: 'molar' or 'mass' (optional, 'molar' as standard)
            specify whether molar or mass fractions are specified.
        """

        Material.__init__(self)

        assert(len(phases) > 0)
        self.phases = phases

        if fractions is not None:
            self.set_fractions(fractions, fraction_type)
        else:
            self.molar_fractions = None

        self.set_averaging_scheme('VoigtReussHill')
        self.name=name
        
    def __str__(self):
        string='Composite: {0}'.format(self.name)
        try:
            string += '\n  P, T: {0:.4g} Pa, {1:.4g} K'.format(self.pressure, self.temperature)
        except:
            pass
        string+='\nPhase and endmember fractions:'
        for phase, fraction in zip(*self.unroll()):
            string+='\n  {0}: {1}'.format(phase.name, fraction)
            if isinstance(phase, SolidSolution):
               for i in range(phase.n_endmembers): 
                   string+='\n    {0}: {1}'.format(phase.endmember_names[i], phase.molar_fractions[i])
        return string

    def set_fractions(self, fractions, fraction_type='molar'):
        """
        Change the fractions of the phases of this Composite.
        Resets cached properties

        Parameters
        ----------
        fractions: list of floats
            molar or mass fraction for each phase.
        fraction_type: 'molar' or 'mass'
            specify whether molar or mass fractions are specified.
        """
        assert(len(self.phases) == len(fractions))

        try:
            total = sum(fractions)
        except TypeError:
            raise Exception(
                "Since v0.8, burnman.Composite takes an array of Materials, then an array of fractions")

        for f in fractions:
            assert (f >= -1e-12)

        self.reset()

        if abs(total - 1.0) > 1e-12:
            warnings.warn(
                "Warning: list of fractions does not add up to one but %g. Normalizing." % total)
            corrected_fractions = [fr / total for fr in fractions]
            fractions = corrected_fractions

        if fraction_type == 'molar':
            molar_fractions = fractions
        elif fraction_type == 'mass':
            molar_fractions = self._mass_to_molar_fractions(
                self.phases, fractions)
        else:
            raise Exception(
                "Fraction type not recognised. Please use 'molar' or mass")

        # Set minimum value of a molar fraction at 0.0 (rather than -1.e-12)
        self.molar_fractions = [max(0.0, fraction)
                                for fraction in molar_fractions]

    def set_method(self, method):
        """
        set the same equation of state method for all the phases in the composite
        """
        for phase in self.phases:
            phase.set_method(method)
        # Clear the cache on resetting method
        self.reset()

    def set_averaging_scheme(self, averaging_scheme):
        """
        Set the averaging scheme for the moduli in the composite.
        Default is set to VoigtReussHill, when Composite is initialized.
        """

        if type(averaging_scheme) == str:
            self.averaging_scheme = getattr(
                averaging_schemes, averaging_scheme)()
        else:
            self.averaging_scheme = averaging_scheme
        # Clear the cache on resetting averaging scheme
        self.reset()

    def set_state(self, pressure, temperature):
        """
        Update the material to the given pressure [Pa] and temperature [K].
        """
        Material.set_state(self, pressure, temperature)
        for phase in self.phases:
            phase.set_state(pressure, temperature)

    def debug_print(self, indent=""):
        print("{0}Composite: {1}".format(indent, self.name))
        indent += "  "
        if self.molar_fractions is None:
            for i, phase in enumerate(self.phases):
                phase.debug_print(indent + "  ")
        else:
            for i, phase in enumerate(self.phases):
                print("%s%g of" % (indent, self.molar_fractions[i]))
                phase.debug_print(indent + "  ")

    def unroll(self):
        if self.molar_fractions is None:
            raise Exception(
                "Unroll only works if the composite has defined fractions.")
        phases = []
        fractions = []
        for i, phase in enumerate(self.phases):
            p_mineral, p_fraction = phase.unroll()
            check_pairs(p_mineral, p_fraction)
            fractions.extend([f * self.molar_fractions[i] for f in p_fraction])
            phases.extend(p_mineral)
        return phases, fractions

    def to_string(self):
        """
        return the name of the composite
        """
        return "{0}: {1}".format(self.__class__.__name__, self.name)

    @material_property
    def element_set(self):
        """
        Returns a list of elements which could be contained in the composite.
        This composite must be composed only of instances of the Mineral
        or SolidSolution class.
        """
        
        elements = []
        for i, ph in enumerate(self.phases):
            if isinstance(ph, SolidSolution):
                for f in ph.endmember_formulae:
                    elements.extend(f.keys())
            elif isinstance(ph, Mineral):
                elements.extend(ph.formula.keys())
            else:
                raise Exception('Composite includes a phase which is neither a mineral nor a solid solution.')

        return list(set(elements))


    @material_property
    def endmembers_per_phase(self):
        """
        Returns a list of integers corresponding to the number of
        endmembers per phase in the Composite (1 for a Mineral).
        This composite must be composed only of instances of the
        Mineral or SolidSolution class.
        """
        n_mbrs = []
        for ph in self.phases:
            if isinstance(ph, SolidSolution):
                n_mbrs.append(len(ph.endmembers))
            elif isinstance(ph, Mineral):
                n_mbrs.append(1)
            else:
                raise Exception('Composite includes a phase which is neither an instance of a burnman.Mineral nor a burnman.SolidSolution.')
        return n_mbrs
        
    @material_property
    def compositional_constraints(self):
        """
        Returns a matrix of the compositional limits (non-negativity constraints)
        for the composite material. Each row constitutes a non-negativity constraint.
        This composite must be composed only of instances of the
        Mineral or SolidSolution class.
        """
        site_occupancies = np.empty([0, 0])
        for ph in self.phases:
            n_pad = site_occupancies.shape[1]

            if isinstance(ph, SolidSolution): 
                if site_occupancies.shape[1] == 0:
                    site_occupancies = np.empty([0, ph.solution_model.endmember_occupancies.shape[1]])
                else:
                    site_occupancies = np.pad(site_occupancies,
                                              ((0, 0),
                                               (0, ph.solution_model.endmember_occupancies.shape[1])),
                                              'constant', constant_values=((0, 0), (0, 0)))
                
                site_occupancies = np.vstack((site_occupancies,
                                              np.pad(ph.solution_model.endmember_occupancies,
                                                     ((0, 0), (n_pad,0)),
                                                     'constant', constant_values=(0., 0.))))
        
            elif isinstance(ph, Mineral):
                site_occupancies = np.pad(site_occupancies,
                                          ((0, 1), (0, 1)),
                                          'constant',
                                          constant_values=((0, 0), (0, 0)))
                site_occupancies[-1,-1] = 1.
            else:
                raise Exception('Composite includes a phase which is neither an instance of a burnman.Mineral nor a burnman.SolidSolution.')
        
        return site_occupancies.T
            
    def stoichiometric_matrix(self, elements = None, excluded_endmembers = None,
                              calculate_subspaces = False):
        """
        Create a 2D matrix containing the number of atoms of element j in
        endmember i of all phases in the composite. Elements to be included and
        endmembers to be excluded are passed as optional arguments to the function.

        Optionally, this function also returns a sparse basis for
        the left null space and column space of the stoichiometric matrix.

        This composite must be composed only of instances of the Mineral
        or SolidSolution class.

        Parameters
        ----------
        elements : list of strings
            list of elements from which to construct the stoichiometric matrix.
            (default = None, in which case the element list is calculated from
            the function element_set)
        excluded_endmembers : list of lists of integers
            indices of the endmembers in each phase to be excluded
            from the stoichiometric matrix. For example, in a composite
            made of four phases, where the last is a Mineral instance and the
            others are SolidSolution instances, the second endmember of the
            first phase and the last phase can be excluded thus:
            excluded_endmembers = [[1], [], [], [0]]
            (default = None)
        calculate_subspaces : bool
            If true, function returns the stoichiometric matrix
            and its left nullspace and columnspace as a tuple

        Returns
        -------
        stoichiometric_matrix : 2D numpy array
            The stoichiometric matrix. Returned either
            as the only object (if calculate_subspaces = False),
            or the first within a tuple of length 3
            (if calculate_subspaces = True).
        S_leftnullspace : 2D numpy array
            The left null space of the stoichiometric matrix
            (returned only if calculate_subspaces = True)
        S_columnspace : 2D numpy array
            The column space of the stoichiometric matrix
            (returned only if calculate_subspaces = True)
        """
        
        def matrix_func(i,j):
            e = elements[j]
            if e in formulae[i]:
                return nsimplify(formulae[i][e])
            else:
                return 0

        if elements is None:
            elements = self.element_set
                    
        if excluded_endmembers is None:
            excluded_endmembers = [[]]*len(self.phases)

        formulae = []
        for i, ph in enumerate(self.phases):
            if isinstance(ph, SolidSolution):
                formulae.extend([f for j, f in enumerate(ph.endmember_formulae)
                                 if j not in excluded_endmembers[i]])
            elif isinstance(ph, Mineral):
                if excluded_endmembers is not [0]:
                    formulae.append(ph.formula)
            else:
                raise Exception('Composite includes a phase which is neither an instance of a burnman.Mineral nor a burnman.SolidSolution.')
        
        stoichiometric_matrix = Matrix( len(formulae), len(elements), matrix_func )

        if calculate_subspaces:
            S_leftnullspace = np.array(stoichiometric_matrix.T.nullspace()).astype(float)
            S_columnspace = np.array(stoichiometric_matrix.columnspace()).astype(float)
            stoichiometric_matrix = np.array(stoichiometric_matrix).astype(float)
            return (stoichiometric_matrix, S_leftnullspace, S_columnspace)
        else:
            stoichiometric_matrix = np.array(stoichiometric_matrix).astype(float)
            return stoichiometric_matrix

    def set_potential_phase_amounts_from_bulk(self, elemental_composition):
        """
        Sets the phase proportions and total number of moles
        of the composite to minimize (using non-negative least squares)
        the difference between the composite composition and 
        elemental composition. 
        Caution: The solution will not necessarily be unique.

        Parameters
        ----------
        elemental_composition : dictionary
            A dictionary containing the number of atoms of each element
            in the composition

        Returns
        -------
        amounts : numpy array
            The inverted amounts of each phase in the assemblage
        residuals : dictionary
            A dictionary of the residuals of all the element amounts
        rnorm : float
            The overall residual || residuals ||_2.
        """

        formulae = [phase.formula for phase in self.phases]
        elements = list(set().union(*formulae+[elemental_composition]))
        
        b = np.array([elemental_composition[e]
                      if e in elemental_composition else 0.
                      for e in elements])
        
        A = np.array([[f[e]
                       if e in f else 0.
                       for f in formulae]
                      for e in elements])
        
        x, rnorm = opt.nnls(A, b)
        res_array = b - A.dot(x)
        residuals = {e: res_array[i] for i, e in enumerate(elements)}

        self.n_moles = sum(x)
        self.set_fractions(x/self.n_moles)
        return (x, residuals, rnorm)
        
        
    def set_potential_composition_from_bulk(self, elemental_composition,
                                            unfitted_elements=None,
                                            excluded_endmembers=None,
                                            exact_solution_required = True,
                                            use_solution_guesses = True):
        """
        Sets phase fractions and compositions of solutions within the
        composite to values which minimize (in a least-squares sense) the
        difference between a provided bulk composition and the composition
        of the assemblage. Also adds/modifies the composite attribute 
        n_moles, which corresponds to the number of moles of the composite
        required to convert the phase fractions into phase amounts.

        A number of parameters provide flexibility in the minimization
        procedure. All elements not in elemental_composition or
        unfitted_elements will be set to zero. All endmembers not in
        excluded_endmembers will be considered.

        In addition to setting the phase fractions and compositions of the
        composite, the function also returns the endmember proportions and
        residuals as numpy arrays.

        Parameters
        ----------
        elemental_composition : dictionary
            A dictionary containing the number of atoms of each element
            in the composition
        unfitted_elements : list of strings
            A list of the elements not to be included in the fitting.
            For example, oxygen can be ignored if the oxidation state
            is unknown.
        excluded_endmembers : list of lists of integers
            indices of the endmembers in each phase to be excluded
            from the stoichiometric matrix. For example, in a composite
            made of four phases, where the last is a Mineral instance and the
            others are SolidSolution instances, the second endmember of the
            first phase and the last phase can be excluded thus:
            excluded_endmembers = [[1], [], [], [0]]
            (default = None)
        exact_solution_required : bool
            If true, the composition is applied as a set of equality
            constraints.
            (default = True)
        use_solution_guesses : bool
            If true, the solutions are checked for the guess attribute.
            If the attribute is present, the function attempts to minimize
            the euclidean distance between the guessed composition(s) and
            the optimized composition(s). If true, exact_solution_required
            must also be True.
            (default = True)
        Returns
        -------
        endmember_amounts : numpy array
            A numpy array of the amounts of each endmember 
        residuals : dictionary
            A dictionary containing the residuals in each of the elements.
        rnorm : float
            The overall residual || residuals ||_2.
        """
        
        if use_solution_guesses and not exact_solution_required:
            raise Exception('If use_solution_guesses is True, '
                            'exact_solution_required must also be set to True.')
            
        if excluded_endmembers is None:
            excluded_endmembers = [[]]*len(self.phases)

        elements = self.element_set
        if unfitted_elements is not None:
            elements = [e for e in elements if e not in unfitted_elements]

        bulk_composition_vector = np.array([elemental_composition[e]
                                            if e in elemental_composition
                                            else 0
                                            for e in elements])

        excluded_indices = []
        m = 0
        for i, n in enumerate(self.endmembers_per_phase):
            excluded_indices.extend([m + v for v in excluded_endmembers[i]])
            m += n
        
        S = self.stoichiometric_matrix(elements, excluded_endmembers)
        E = np.delete(self.compositional_constraints.T, excluded_indices, axis=0)
        Enull = np.array(Matrix(E).nullspace()).astype(float)
        
        A = S.T.dot(pinv(E.T))
    
        # First, let's make a guess at a composition using nnls
        # This will not satisfy the site occupancy constraints,
        # but should provide a reasonable starting guess for the constrained minimization
        Aprime = np.vstack((A, Enull))
        bprime = np.concatenate((bulk_composition_vector, np.zeros(len(Enull))))
        xprime_guess, residual = opt.nnls(Aprime, bprime)

        # Now, let's run a constrained minimization
        # Define constraints and bounds
        # We need some python magic to do this dynamically
        endmember_constraints = lambda nullspace: [{'type': 'eq',
                                                    'fun': lambda x, eq=eq: eq.dot(x)}
                                                   for eq in nullspace]

        cons = endmember_constraints(Enull)

        if exact_solution_required:
            bulk_constraints = lambda A, b: [{'type': 'eq',
                                              'fun': lambda x, Ai=A[i], bi=b[i]: Ai.dot(x) - bi}
                                             for i in range(len(A))]
            cons.extend(bulk_constraints(A, bulk_composition_vector))

        bounds = [[0., None]]*len(E.T)

        fn = lambda x, A, b: np.linalg.norm(A.dot(x) - b)
        sol = opt.minimize(fn, xprime_guess, args=(A, bulk_composition_vector),
                           method='SLSQP', bounds=bounds, constraints=cons)


        # Optionally, use this solution and refine based on some
        # solution starting guesses
        if use_solution_guesses:
            guesses = []
            
            n_total = sum(self.endmembers_per_phase)
            mul = np.zeros((n_total, n_total))

            m = 0
            guesses_exist = False
            for i, n in enumerate(self.endmembers_per_phase):
                if hasattr(self.phases[i], 'guess'):
                    guesses_exist = True
                    guesses.extend(list(self.phases[i].guess))
                else:
                    guesses.extend([0.]*n)
                    
                mul[m:m+n, m:m+n] = 1.
                m += n

            if guesses_exist:
                guesses = np.delete(np.array(guesses), excluded_indices, axis=0)
                mul = np.delete(mul, excluded_indices, axis=0)
                mul = np.delete(mul, excluded_indices, axis=1)
                
                def fn(x, E, mul):
                    mbr_amounts = pinv(E.T).dot(sol.x)
                    phase_amounts = np.array([max(1.e-10, v) for v in mul.dot(mbr_amounts)])

                    mbr_fractions = mbr_amounts/phase_amounts
                    guessed_fractions = np.copy(mbr_fractions)

                    for i, phase_amount in enumerate(mul.dot(guesses)):
                        if phase_amount > 1.e-10:
                            guessed_fractions[i] = guesses[i]/phase_amount

                    return np.linalg.norm(mbr_fractions - guessed_fractions)

                solg = opt.minimize(fn, sol.x, args=(E, mul),
                                    method='SLSQP', bounds=bounds, constraints=cons)

                if solg.success:
                    sol = solg
                else:
                    warnings.warn('Solver failed to minimize the residuals in composition. '
                                  'Falling back to solution ignoring compositional guesses.')
                        
            else:
                warnings.warn('No compositional guesses exist for this Composite. Set the function variable use_solution_guesses = False to remove this warning.', stacklevel=2)

        endmember_amounts = pinv(E.T).dot(sol.x)

        # Now calcuate residuals and insert solution into the composite
        if sol.success or sol.fun < 1.e-12:
            
            residual_values = S.T.dot(endmember_amounts) - bulk_composition_vector
            residuals = {e: residual_values[i] for i, e in enumerate(elements)}
            rnorm = np.sqrt(sum(residual_values*residual_values))
            
            for idx in excluded_indices:
                endmember_amounts = np.insert(endmember_amounts, idx, 0.)

            m = 0
            n_moles = []
            for i, n in enumerate(self.endmembers_per_phase):
                if n > 1:
                    ci = endmember_amounts[m:m+n]
                    phase_amount = sum(ci)
                    endmember_proportions = ci/sum(ci)

                    if phase_amount < 1.e-10:
                        phase_amount = 0.
                        if hasattr(self.phases[i], 'guess'):
                            endmember_proportions = self.phases[i].guess
                        else:
                            endmember_proportions = [0]*len(endmember_proportions)
                            endmember_proportions[0] = 1.

                    else:
                        # check validity of composition. If not completely valid, run a last minimization
                        constraints = self.compositional_constraints[:,m:m+n].dot(endmember_proportions)
                        if not all(constraints >= 0.):
                            endmember_constraints = lambda constraints: [{'type': 'ineq',
                                                                          'fun': lambda x, c=c: 1.e20*c.dot(x)}
                                                                         for c in constraints]
                            
                            cons = endmember_constraints(self.compositional_constraints[:,m:m+n])
        
                            fn = lambda x, x0: np.linalg.norm(x - x0)
                            solv = opt.minimize(fn, endmember_proportions, args=(endmember_proportions),
                                                  method='SLSQP', constraints=cons)
                            endmember_proportions = solv.x
                            
                        endmember_proportions[np.abs(endmember_proportions) < 1.e-12] = 0.
                    n_moles.append(phase_amount)
                    self.phases[i].set_composition(endmember_proportions)
                else:
                    n_moles.append(endmember_amounts[m])
                m += n
            self.n_moles = sum(n_moles)
            self.set_fractions(np.array(n_moles)/self.n_moles)
                    
            return (endmember_amounts, residuals, rnorm)
        else:
            raise Exception("A solution was not found.")
        

    @material_property
    def molar_internal_energy(self):
        """
        Returns molar internal energy of the mineral [J/mol]
        Aliased with self.energy
        """
        U = sum(phase.molar_internal_energy * molar_fraction for (
                phase, molar_fraction) in zip(self.phases, self.molar_fractions))
        return U

    @material_property
    def molar_gibbs(self):
        """
        Returns molar Gibbs free energy of the composite [J/mol]
        Aliased with self.gibbs
        """
        G = sum(phase.molar_gibbs * molar_fraction for (phase, molar_fraction)
                in zip(self.phases, self.molar_fractions))
        return G

    @material_property
    def molar_helmholtz(self):
        """
        Returns molar Helmholtz free energy of the mineral [J/mol]
        Aliased with self.helmholtz
        """
        F = sum(phase.molar_helmholtz * molar_fraction for (
                phase, molar_fraction) in zip(self.phases, self.molar_fractions))
        return F

    @material_property
    def molar_volume(self):
        """
        Returns molar volume of the composite [m^3/mol]
        Aliased with self.V
        """
        volumes = np.array(
            [phase.molar_volume * molar_fraction for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)])
        return np.sum(volumes)

    @material_property
    def molar_mass(self):
        """
        Returns molar mass of the composite [kg/mol]
        """
        return sum([phase.molar_mass * molar_fraction for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)])

    @material_property
    def density(self):
        """
        Compute the density of the composite based on the molar volumes and masses
        Aliased with self.rho
        """
        densities = np.array([phase.density for phase in self.phases])
        volumes = np.array(
            [phase.molar_volume * molar_fraction for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)])
        return self.averaging_scheme.average_density(volumes, densities)

    @material_property
    def molar_entropy(self):
        """
        Returns enthalpy of the mineral [J]
        Aliased with self.S
        """
        S = sum(phase.molar_entropy * molar_fraction for (
                phase, molar_fraction) in zip(self.phases, self.molar_fractions))
        return S

    @material_property
    def molar_enthalpy(self):
        """
        Returns enthalpy of the mineral [J]
        Aliased with self.H
        """
        H = sum(phase.molar_enthalpy * molar_fraction for (
                phase, molar_fraction) in zip(self.phases, self.molar_fractions))
        return H

    @material_property
    def isothermal_bulk_modulus(self):
        """
        Returns isothermal bulk modulus of the composite [Pa]
        Aliased with self.K_T
        """
        V_frac = np.array([phase.molar_volume * molar_fraction for (
                           phase, molar_fraction) in zip(self.phases, self.molar_fractions)])
        K_ph = np.array(
            [phase.isothermal_bulk_modulus for phase in self.phases])
        G_ph = np.array([phase.shear_modulus for phase in self.phases])

        return self.averaging_scheme.average_bulk_moduli(V_frac, K_ph, G_ph)

    @material_property
    def adiabatic_bulk_modulus(self):
        """
        Returns adiabatic bulk modulus of the mineral [Pa]
        Aliased with self.K_S
        """
        V_frac = np.array([phase.molar_volume * molar_fraction for (
                           phase, molar_fraction) in zip(self.phases, self.molar_fractions)])
        K_ph = np.array(
            [phase.adiabatic_bulk_modulus for phase in self.phases])
        G_ph = np.array([phase.shear_modulus for phase in self.phases])

        return self.averaging_scheme.average_bulk_moduli(V_frac, K_ph, G_ph)

    @material_property
    def isothermal_compressibility(self):
        """
        Returns isothermal compressibility of the composite (or inverse isothermal bulk modulus) [1/Pa]
        Aliased with self.beta_T
        """
        return 1. / self.isothermal_bulk_modulus

    @material_property
    def adiabatic_compressibility(self):
        """
        Returns isothermal compressibility of the composite (or inverse isothermal bulk modulus) [1/Pa]
        Aliased with self.beta_S
        """
        return 1. / self.adiabatic_bulk_modulus

    @material_property
    def shear_modulus(self):
        """
        Returns shear modulus of the mineral [Pa]
        Aliased with self.G
        """
        V_frac = np.array([phase.molar_volume * molar_fraction for (
                           phase, molar_fraction) in zip(self.phases, self.molar_fractions)])
        K_ph = np.array(
            [phase.adiabatic_bulk_modulus for phase in self.phases])
        G_ph = np.array([phase.shear_modulus for phase in self.phases])

        return self.averaging_scheme.average_shear_moduli(V_frac, K_ph, G_ph)

    @material_property
    def p_wave_velocity(self):
        """
        Returns P wave speed of the composite [m/s]
        Aliased with self.v_p
        """
        return np.sqrt((self.adiabatic_bulk_modulus + 4. / 3. *
                        self.shear_modulus) / self.density)

    @material_property
    def bulk_sound_velocity(self):
        """
        Returns bulk sound speed of the composite [m/s]
        Aliased with self.v_phi
        """
        return np.sqrt(self.adiabatic_bulk_modulus / self.density)

    @material_property
    def shear_wave_velocity(self):
        """
        Returns shear wave speed of the composite [m/s]
        Aliased with self.v_s
        """
        return np.sqrt(self.shear_modulus / self.density)

    @material_property
    def grueneisen_parameter(self):
        """
        Returns grueneisen parameter of the composite [unitless]
        Aliased with self.gr
        """
        return self.thermal_expansivity * self.isothermal_bulk_modulus * self.molar_volume / self.molar_heat_capacity_v

    @material_property
    def thermal_expansivity(self):
        """
        Returns thermal expansion coefficient of the composite [1/K]
        Aliased with self.alpha
        """
        volumes = np.array(
            [phase.molar_volume * molar_fraction for (phase, molar_fraction) in zip(self.phases, self.molar_fractions)])
        alphas = np.array([phase.thermal_expansivity for phase in self.phases])
        return self.averaging_scheme.average_thermal_expansivity(volumes, alphas)

    @material_property
    def molar_heat_capacity_v(self):
        """
        Returns molar_heat capacity at constant volume of the composite [J/K/mol]
        Aliased with self.C_v
        """
        c_v = np.array([phase.molar_heat_capacity_v for phase in self.phases])
        return self.averaging_scheme.average_heat_capacity_v(self.molar_fractions, c_v)

    @material_property
    def molar_heat_capacity_p(self):
        """
        Returns molar heat capacity at constant pressure of the composite [J/K/mol]
        Aliased with self.C_p
        """
        c_p = np.array([phase.molar_heat_capacity_p for phase in self.phases])
        return self.averaging_scheme.average_heat_capacity_p(self.molar_fractions, c_p)

    def _mass_to_molar_fractions(self, phases, mass_fractions):
        """
        Converts a set of mass fractions for phases into a set of molar fractions.

        Parameters
        ----------
        phases : list of :class:`burnman.Material`
        The list of phases for which fractions should be converted.

        mass_fractions : list of floats
        The list of mass fractions of the input phases.

        Returns
        -------
        molar_fractions : list of floats
        The list of molar fractions corresponding to the input molar fractions
        """
        total_moles = sum(
            mass_fraction / phase.molar_mass for mass_fraction, phase in zip(mass_fractions, phases))
        molar_fractions = [mass_fraction / (phase.molar_mass * total_moles)
                           for mass_fraction, phase in zip(mass_fractions, phases)]
        return molar_fractions
