import unittest
import os, sys
sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals


# TODO: test composite that changes number of entries

class test(unittest.TestCase):
    
    def test_simple(self):
        inp1 = {'Mg':0.213, 'Fe': 0.0626, 'Si':0.242, 'Ca':0., 'Al':0.} # wt%
        phase_per,rel_mol_per = burnman.calculate_phase_percents(inp1) 
        
        StartP = 23.83 #in GPa
        EndP = 110.0
        deltaP = 1.

        #P,T,a,b,frac_mol_pv,frac_mol_mw    = part_coef_calc(inp2,StartP,EndP,deltaP)
 
        gt = lambda p: burnman.geotherm.brown_shankland(p)
        pressure = StartP
        temperature = gt([StartP])[0]
        (a,b) = burnman.calculate_partition_coefficient(pressure, temperature, rel_mol_per, 0.5)
        self.assertAlmostEqual(a, 0.14347274523)
        self.assertAlmostEqual(b, 0.108980503625)
        


if __name__ == '__main__':
    unittest.main()
