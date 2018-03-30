from __future__ import absolute_import
from __future__ import print_function
import unittest
import os
import sys

sys.path.insert(1, os.path.abspath('..'))
import warnings

import burnman
from burnman import minerals
from burnman.processchemistry import dictionarize_formula, formula_mass

from util import BurnManTest


class mypericlase(burnman.Mineral):

    """
    Stixrude & Lithgow-Bertelloni 2005 and references therein
    """

    def __init__(self):
        self.params = {
            'equation_of_state': 'slb3',
            'T_0': 300.,
            'P_0': 0.,
            'V_0': 11.24e-6,
            'K_0': 161.0e9,
            'Kprime_0': 3.8,
            'G_0': 131.0e9,
            'Gprime_0': 2.1,
            'molar_mass': .0403,
            'n': 2,
            'Debye_0': 773.,
            'grueneisen_0': 1.5,
            'q_0': 1.5,
            'eta_s_0': 2.8}


class Fe_Dewaele(burnman.Mineral):

    """
    Dewaele et al., 2006, Physical Review Letters
    """

    def __init__(self):
        self.params = {
            'equation_of_state': 'vinet',
            'P_0': 0.,
            'V_0': 6.75e-6,
            'K_0': 163.4e9,
            'Kprime_0': 5.38,
            'molar_mass': 0.055845,
            'n': 1}


class Liquid_Fe_Anderson(burnman.Mineral):

    """
    Anderson & Ahrens, 1994 JGR
    """

    def __init__(self):
        self.params = {
            'equation_of_state': 'bm4',
            'P_0': 0.,
            'V_0': 7.95626e-6,
            'K_0': 109.7e9,
            'Kprime_0': 4.66,
            'Kprime_prime_0': -0.043e-9,
            'molar_mass': 0.055845,
        }


class outer_core_rkprime(burnman.Mineral):

    """
    Stacey and Davis, 2004 PEPI (Table 5)
    """

    def __init__(self):
        self.params = {
            'equation_of_state': 'rkprime',
            'P_0': 0.,
            'V_0': 0.055845/6562.54,
            'K_0': 124.553e9,
            'Kprime_0': 4.9599,
            'Kprime_inf': 3.0,
            'molar_mass': 0.055845,
        }

        
class periclase_morse(burnman.Mineral):

    """
    Periclase parameters from SLB dataset (which uses BM3)
    """

    def __init__(self):
        formula = 'MgO'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'Periclase',
            'formula': formula,
            'equation_of_state': 'morse',
            'P_0': 0.,
            'V_0': 1.1244e-05,
            'K_0': 1.613836e+11,
            'Kprime_0': 3.84045,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula)}
        
        

class eos(BurnManTest):

    def test_reference_values(self):
        rock = mypericlase()
        pressure = 0.
        temperature = 300.
        eoses = [
            burnman.eos.SLB2(), burnman.eos.SLB3(), burnman.eos.BM2(), burnman.eos.BM3()]

        for i in eoses:
            Volume_test = i.volume(pressure, temperature, rock.params)
            self.assertFloatEqual(Volume_test, rock.params['V_0'])
            Kt_test = i.isothermal_bulk_modulus(
                pressure, 300., rock.params['V_0'], rock.params)
            self.assertFloatEqual(Kt_test, rock.params['K_0'])
            # K_S is based on 0 reference temperature:
            Kt_test = i.isothermal_bulk_modulus(
                pressure, 0., rock.params['V_0'], rock.params)
            K_test = i.adiabatic_bulk_modulus(
                pressure, 0., rock.params['V_0'], rock.params)
            self.assertFloatEqual(K_test, Kt_test)
            G_test = i.shear_modulus(
                pressure, temperature, rock.params['V_0'], rock.params)
            self.assertFloatEqual(G_test, rock.params['G_0'])
            Density_test = i.density(rock.params['V_0'], rock.params)
            self.assertFloatEqual(
                Density_test, rock.params['molar_mass'] / rock.params['V_0'])
            alpha_test = i.thermal_expansivity(
                pressure, temperature, rock.params['V_0'], rock.params)
            Cp_test = i.molar_heat_capacity_p(
                pressure, temperature, rock.params['V_0'], rock.params)
            Cv_test = i.molar_heat_capacity_v(
                pressure, temperature, rock.params['V_0'], rock.params)
            Grun_test = i.grueneisen_parameter(
                pressure, temperature, rock.params['V_0'], rock.params)

        eoses_thermal = [burnman.eos.SLB2(), burnman.eos.SLB3()]
        for i in eoses_thermal:
            Cp_test = i.molar_heat_capacity_p(
                pressure, temperature, rock.params['V_0'], rock.params)
            self.assertFloatEqual(Cp_test, 37.076768469502042)
            Cv_test = i.molar_heat_capacity_v(
                pressure, temperature, rock.params['V_0'], rock.params)
            self.assertFloatEqual(Cv_test, 36.577717628901553)
            alpha_test = i.thermal_expansivity(
                pressure, temperature, rock.params['V_0'], rock.params)
            self.assertFloatEqual(alpha_test, 3.031905596878513e-05)
            Grun_test = i.grueneisen_parameter(
                pressure, temperature, rock.params['V_0'], rock.params)
            self.assertFloatEqual(Grun_test, rock.params['grueneisen_0'])

    def test_reference_values_vinet(self):
        rock = Fe_Dewaele()
        pressure = 0.
        temperature = 300.
        eos = burnman.eos.Vinet()

        Volume_test = eos.volume(pressure, temperature, rock.params)
        self.assertFloatEqual(Volume_test, rock.params['V_0'])
        Kt_test = eos.isothermal_bulk_modulus(
            pressure, 300., rock.params['V_0'], rock.params)
        self.assertFloatEqual(Kt_test, rock.params['K_0'])
        Density_test = eos.density(rock.params['V_0'], rock.params)
        self.assertFloatEqual(
            Density_test, rock.params['molar_mass'] / rock.params['V_0'])

    def test_reference_values_bm4(self):
        rock = Liquid_Fe_Anderson()
        pressure = 0.
        temperature = 300.
        eos = burnman.eos.BM4()

        Volume_test = eos.volume(pressure, temperature, rock.params)
        self.assertFloatEqual(Volume_test, rock.params['V_0'])
        Kt_test = eos.isothermal_bulk_modulus(
            pressure, 300., rock.params['V_0'], rock.params)
        self.assertFloatEqual(Kt_test, rock.params['K_0'])
        Density_test = eos.density(rock.params['V_0'], rock.params)
        self.assertFloatEqual(
            Density_test, rock.params['molar_mass'] / rock.params['V_0'])

    def test_reference_values_morse(self):
        rock = periclase_morse()
        pressure = 0.
        temperature = 300.
        eos = burnman.eos.Morse()
        Volume_test = eos.volume(pressure, temperature, rock.params)
        self.assertFloatEqual(Volume_test, rock.params['V_0'])
        Kt_test = eos.isothermal_bulk_modulus(
            pressure, 300., rock.params['V_0'], rock.params)
        self.assertFloatEqual(Kt_test, rock.params['K_0'])
        Density_test = eos.density(rock.params['V_0'], rock.params)
        self.assertFloatEqual(
            Density_test, rock.params['molar_mass'] / rock.params['V_0'])
        
    def test_reference_values_rkprime(self):
        rock = outer_core_rkprime()
        pressure = 0.
        temperature = 300.
        eos = burnman.eos.RKprime()
        Volume_test = eos.volume(pressure, temperature, rock.params)
        self.assertFloatEqual(Volume_test, rock.params['V_0'])
        Kt_test = eos.isothermal_bulk_modulus(
            pressure, 300., rock.params['V_0'], rock.params)
        self.assertFloatEqual(Kt_test, rock.params['K_0'])
        Density_test = eos.density(rock.params['V_0'], rock.params)
        self.assertFloatEqual(
            Density_test, rock.params['molar_mass'] / rock.params['V_0'])

    def test_reference_values_aa(self):
        rock = minerals.other.liquid_iron()
        pressure = rock.params['P_0']
        temperature = rock.params['T_0']
        eos = burnman.eos.AA()
        Volume_test = eos.volume(pressure, temperature, rock.params)
        self.assertFloatEqual(Volume_test, rock.params['V_0'])
        Ks_test = eos.adiabatic_bulk_modulus(
            pressure, temperature, rock.params['V_0'], rock.params)
        self.assertFloatEqual(Ks_test, rock.params['K_S'])
        Density_test = eos.density(rock.params['V_0'], rock.params)
        self.assertFloatEqual(
            Density_test, rock.params['molar_mass'] / rock.params['V_0'])
        
        
class test_eos_validation(BurnManTest):

    def test_no_shear_error(self):
        # The validation should place nans in for the shear parameters
        # If any exceptions or warnings are raised, fail.
        class mymineralwithoutshear(burnman.Mineral):

            def __init__(self):
                self.params = {
                    'equation_of_state': 'slb3',
                    'V_0': 11.24e-6,
                    'K_0': 161.0e9,
                    'Kprime_0': 3.8,
                    'molar_mass': .0403,
                    'n': 2,
                    'Debye_0': 773.,
                    'grueneisen_0': 1.5,
                    'q_0': 1.5,
                    'eta_s_0': 2.8}
                burnman.Mineral.__init__(self)

        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            # Trigger warning
            shearless = mymineralwithoutshear()
            if len(w) != 0:
                self.fail("Caught unexpected warning: " + str(w[-1]))
            try:
                x = shearless.params['G_0']
                y = shearless.params['Gprime_0']
                z = shearless.params['F_0']
                z = shearless.params['eta_s_0']
            except KeyError:
                self.fail('Parameter padding failed in validation')
                pass

    def test_dumb_parameter_values(self):

        class mymineralwithnegativekprime(burnman.Mineral):

            def __init__(self):
                self.params = {
                    'equation_of_state': 'slb3',
                    'V_0': 11.24e-6,
                    'K_0': 161.0e9,
                    'Kprime_0': -4.,
                    'molar_mass': .0403,
                    'n': 2,
                    'Debye_0': 773.,
                    'grueneisen_0': 1.5,
                    'q_0': 1.5,
                    'eta_s_0': 2.8}
                burnman.Mineral.__init__(self)

        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            # Trigger warning
            negative_Kprime = mymineralwithnegativekprime()
            if len(w) == 0:
                print(negative_Kprime.params)
                self.fail(
                    "Did not catch expected warning for negative K prime")

        class mymineralwithkingigapascals(burnman.Mineral):

            def __init__(self):
                self.params = {
                    'equation_of_state': 'slb3',
                    'V_0': 11.24e-6,
                    'K_0': 161.0,
                    'Kprime_0': 3.8,
                    'molar_mass': .0403,
                    'n': 3.14159,
                    'Debye_0': 773.,
                    'grueneisen_0': 1.5,
                    'q_0': 1.5,
                    'eta_s_0': 2.8}
                burnman.Mineral.__init__(self)

        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            # Trigger warning
            low_K = mymineralwithkingigapascals()
            if len(w) == 0:
                self.fail("Did not catch expected warning K in Gpa")

    def test_reference_energies(self):
        m = burnman.Mineral(params={'equation_of_state': 'bm3',
                                    'E_0': 1000.,
                                    'T_0': 100.,
                                    'P_0': 1.e10,
                                    'V_0': 7.95626e-6,
                                    'K_0': 109.7e9,
                                    'Kprime_0': 4.66,
                                    'Kprime_inf': 3.00,
                                    'Kprime_prime_0': -0.043e-9,
                                    'Kdprime_0': -4.66/100.e9,
                                    'molar_mass': 0.055845})

        eoses = ['bm3', 'bm4', 'vinet', 'mt', 'morse', 'rkprime']

        energies = []
        for eos in eoses:
            m.params['equation_of_state'] = eos
            burnman.Mineral.__init__(m)
            m.set_state(m.params['P_0'], m.params['T_0'])
            energies.append(m.molar_internal_energy)
        
        self.assertArraysAlmostEqual(energies, [m.params['E_0']]*len(energies))
    
    def test_energy_derivatives(self):
        m = burnman.Mineral(params={'equation_of_state': 'bm3',
                                    'V_0': 7.95626e-6,
                                    'K_0': 109.7e9,
                                    'Kprime_0': 4.66,
                                    'Kprime_inf': 3.00,
                                    'Kprime_prime_0': -0.043e-9,
                                    'Kdprime_0': -4.66/100.e9,
                                    'molar_mass': 0.055845})

        eoses = ['bm3', 'bm4', 'vinet', 'mt', 'morse', 'rkprime']

        calculated = []
        derivative = []
        for eos in eoses:
            m.params['equation_of_state'] = eos
            burnman.Mineral.__init__(m)

            P_0 = 10.e9
            dP = 1000.
            pressures = [P_0 - 0.5*dP, P_0, P_0 + 0.5*dP]
            temperatures = [0., 0., 0.]

            E, G, H, A, V = m.evaluate(['molar_internal_energy', 'gibbs', 'H', 'helmholtz', 'V'], pressures, temperatures)
            
            calculated.append(P_0)
            derivative.append(-(E[2] - E[0])/(V[2] - V[0]))
            calculated.append(V[1])
            derivative.append((G[2] - G[0])/dP)
            
        self.assertArraysAlmostEqual(calculated, derivative)
                

if __name__ == '__main__':
    unittest.main()
