import unittest
import os, sys
sys.path.insert(1,os.path.abspath('..'))
import warnings

import burnman
from burnman import minerals

from util import BurnManTest


class mypericlase (burnman.Mineral):
    """
    Stixrude & Lithgow-Bertelloni 2005 and references therein 
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
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
            'eta_s_0': 2.8 }

class eos(BurnManTest):

    def test_reference_values(self):
		rock = mypericlase()
		pressure = 0.
		temperature = 300
		#rock.set_method('slb2')
		eoses = [burnman.slb.SLB2(),burnman.slb.SLB3(),burnman.birch_murnaghan.BM2(),burnman.birch_murnaghan.BM3()]
		
		for i in eoses:
            Volume_test = i.volume(pressure, temperature,rock.params)
            self.assertFloatEqual(Volume_test,rock.params['V_0'])
        	K_test = i.adiabatic_bulk_modulus(pressure, temperature,rock.params['V_0'],rock.params)
        	self.assertFloatEqual(K_test,rock.params['K_0'])
        	G_test = i.shear_modulus(pressure, temperature,rock.params['V_0'],rock.params)
        	self.assertFloatEqual(G_test,rock.params['G_0'])
        	Density_test = i.density(pressure, temperature, rock.params)
        	self.assertFloatEqual(Density_test,rock.params['molar_mass']/rock.params['V_0'])
        	Cp_test= i.heat_capacity_p(pressure, temperature,rock.params['V_0'],rock.params)
        	#self.assertFloatEqual(Cp_test,...)
        	Cv_test= i.heat_capacity_v(pressure, temperature,rock.params['V_0'],rock.params)
        	Kt_test = i.isothermal_bulk_modulus(pressure, temperature,rock.params['V_0'],rock.params)
        	self.assertFloatEqual(Kt_test,rock.params['K_0'])
        	alpha_test = i.thermal_expansivity(pressure, temperature,rock.params['V_0'],rock.params)
        	Cp_test= i.heat_capacity_p(pressure, temperature,rock.params['V_0'],rock.params)
        	Cv_test= i.heat_capacity_v(pressure, temperature,rock.params['V_0'],rock.params)
        	Grun_test = i.grueneisen_parameter(pressure, temperature,rock.params['V_0'],rock.params)

		eoses_thermal = [burnman.slb.SLB2(),burnman.slb.SLB3()]
		for i in eoses_thermal:
            Cp_test= i.heat_capacity_p(pressure, temperature,rock.params['V_0'],rock.params)
        	self.assertFloatEqual(Cp_test,37.076768469502042)
        	Cv_test= i.heat_capacity_v(pressure, temperature,rock.params['V_0'],rock.params)
        	self.assertFloatEqual(Cv_test,36.577717628901553)        
        	alpha_test = i.thermal_expansivity(pressure, temperature,rock.params['V_0'],rock.params)
        	self.assertFloatEqual(alpha_test,3.031905596878513e-05)
        	Grun_test = i.grueneisen_parameter(pressure, temperature,rock.params['V_0'],rock.params)
        	self.assertFloatEqual(Grun_test,rock.params['grueneisen_0'])


if __name__ == '__main__':
    unittest.main()
