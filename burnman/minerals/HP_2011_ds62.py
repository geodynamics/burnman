# BurnMan - a lower mantle toolkit
# Copyright (C) 2012-2014, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""

HP_2011_ds62
^^^^^^^^

Minerals from Holland and Powell (2011) and references therein
Update to dataset version 6.2
The values in this document are all in S.I. units, 
unlike those in the original tc-ds62.txt

N.B. VERY IMPORTANT: The excess entropy term in the regular solution model has the opposite sign to the values in Holland and Powell, 2011. This is consistent with its treatment as an excess entropy term (G=H-T*S+P*V), rather than a thermal correction to the interaction parameter (W=W+T*W_T+P*W_P).
"""

from burnman.processchemistry import ProcessChemistry
from burnman.mineral import Mineral
from burnman.solidsolution import SolidSolution

# Mixed Holland and Powell sources (test only)
class garnet(SolidSolution):
    def __init__(self):
        # Name
        self.name='garnet'

        # Endmembers
        base_material = [[py(), '[Mg]3[Al]2Si3O12', 1.0],[alm(), '[Fe]3[Al]2Si3O12', 1.0],[gr(), '[Ca]3[Al]2Si3O12', 2.7],[maj(), '[Mg]3[Mg1/2Si1/2]2Si3O12', 1.0]]

        # Interaction parameters
        excess_enthalpy=[[2.5e3, 29.1e3, 15e3],[10e3,18e3],[48e3]]
        excess_entropy=[[0., 0., 0.],[0., 0.],[0.]]
        excess_volume=[[0., 0.164e-5, 0.],[0., 0.],[0.]]
        
        interaction_parameter=[excess_enthalpy,excess_entropy,excess_volume]
        SolidSolution.__init__(self, base_material, interaction_parameter)



class fo (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg2.0Si1.0O4.0'
       self.params = {
            'name': 'fo',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2172590.0 ,
            'S_0': 95.1 ,
            'V_0': 4.366e-05 ,
            'Cp': [233.3, 0.001494, -603800.0, -1869.7] ,
            'a_0': 2.85e-05 ,
            'K_0': 1.285e+11 ,
            'Kprime_0': 3.84 ,
            'Kdprime_0': -3e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class fa (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe2.0Si1.0O4.0'
       self.params = {
            'name': 'fa',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1477720.0 ,
            'S_0': 151.0 ,
            'V_0': 4.631e-05 ,
            'Cp': [201.1, 0.01733, -1960600.0, -900.9] ,
            'a_0': 2.82e-05 ,
            'K_0': 1.256e+11 ,
            'Kprime_0': 4.68 ,
            'Kdprime_0': -3.7e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class teph (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mn2.0Si1.0O4.0'
       self.params = {
            'name': 'teph',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1733970.0 ,
            'S_0': 155.9 ,
            'V_0': 4.899e-05 ,
            'Cp': [219.6, 0.0, -1292700.0, -1308.3] ,
            'a_0': 2.86e-05 ,
            'K_0': 1.256e+11 ,
            'Kprime_0': 4.68 ,
            'Kdprime_0': -3.7e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class lrn (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca2.0Si1.0O4.0'
       self.params = {
            'name': 'lrn',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2306920.0 ,
            'S_0': 127.6 ,
            'V_0': 5.16e-05 ,
            'Cp': [247.5, -0.003206, 0.0, -2051.9] ,
            'a_0': 2.9e-05 ,
            'K_0': 98500000000.0 ,
            'Kprime_0': 4.07 ,
            'Kdprime_0': -4.1e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'landau_Tc': 1710.0 ,
            'landau_Smax': 10.03 ,
            'landau_Vmax': 5e-07 }

class mont (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca1.0Mg1.0Si1.0O4.0'
       self.params = {
            'name': 'mont',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2251260.0 ,
            'S_0': 109.5 ,
            'V_0': 5.148e-05 ,
            'Cp': [250.7, -0.010433, -797200.0, -1996.1] ,
            'a_0': 2.87e-05 ,
            'K_0': 1.134e+11 ,
            'Kprime_0': 3.87 ,
            'Kdprime_0': -3.4e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class chum (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg9.0Si4.0O18.0H2.0'
       self.params = {
            'name': 'chum',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -9613540.0 ,
            'S_0': 440.5 ,
            'V_0': 0.00019801 ,
            'Cp': [1071.0, -0.016533, -7899600.0, -7373.9] ,
            'a_0': 3.2e-05 ,
            'K_0': 1.199e+11 ,
            'Kprime_0': 4.58 ,
            'Kdprime_0': -3.8e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class chdr (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg5.0Si2.0O10.0H2.0'
       self.params = {
            'name': 'chdr',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -5254890.0 ,
            'S_0': 260.0 ,
            'V_0': 0.00011084 ,
            'Cp': [625.0, -0.001088, -2259900.0, -4910.7] ,
            'a_0': 1.82e-05 ,
            'K_0': 1.161e+11 ,
            'Kprime_0': 4.8 ,
            'Kdprime_0': -4.1e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class mwd (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg2.0Si1.0O4.0'
       self.params = {
            'name': 'mwd',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2138520.0 ,
            'S_0': 93.9 ,
            'V_0': 4.051e-05 ,
            'Cp': [208.7, 0.003942, -1709500.0, -1302.8] ,
            'a_0': 2.37e-05 ,
            'K_0': 1.726e+11 ,
            'Kprime_0': 3.84 ,
            'Kdprime_0': -2.2e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class fwd (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe2.0Si1.0O4.0'
       self.params = {
            'name': 'fwd',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1467900.0 ,
            'S_0': 146.0 ,
            'V_0': 4.321e-05 ,
            'Cp': [201.1, 0.01733, -1960600.0, -900.9] ,
            'a_0': 2.73e-05 ,
            'K_0': 1.69e+11 ,
            'Kprime_0': 4.35 ,
            'Kdprime_0': -2.6e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class mrw (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg2.0Si1.0O4.0'
       self.params = {
            'name': 'mrw',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2127680.0 ,
            'S_0': 90.0 ,
            'V_0': 3.949e-05 ,
            'Cp': [213.3, 0.00269, -1410400.0, -1495.9] ,
            'a_0': 2.01e-05 ,
            'K_0': 1.781e+11 ,
            'Kprime_0': 4.35 ,
            'Kdprime_0': -2.4e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class frw (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe2.0Si1.0O4.0'
       self.params = {
            'name': 'frw',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1471760.0 ,
            'S_0': 140.0 ,
            'V_0': 4.203e-05 ,
            'Cp': [166.8, 0.04261, -1705400.0, -541.4] ,
            'a_0': 2.22e-05 ,
            'K_0': 1.977e+11 ,
            'Kprime_0': 4.92 ,
            'Kdprime_0': -2.5e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class mpv (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg1.0Si1.0O3.0'
       self.params = {
            'name': 'mpv',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1443030.0 ,
            'S_0': 62.6 ,
            'V_0': 2.445e-05 ,
            'Cp': [149.3, 0.002918, -2983000.0, -799.1] ,
            'a_0': 1.87e-05 ,
            'K_0': 2.51e+11 ,
            'Kprime_0': 4.14 ,
            'Kdprime_0': -1.6e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class fpv (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe1.0Si1.0O3.0'
       self.params = {
            'name': 'fpv',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1084640.0 ,
            'S_0': 91.0 ,
            'V_0': 2.548e-05 ,
            'Cp': [133.2, 0.01083, -3661400.0, -314.7] ,
            'a_0': 1.87e-05 ,
            'K_0': 2.81e+11 ,
            'Kprime_0': 4.14 ,
            'Kdprime_0': -1.6e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class apv (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Al2.0O3.0'
       self.params = {
            'name': 'apv',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1646630.0 ,
            'S_0': 51.8 ,
            'V_0': 2.54e-05 ,
            'Cp': [139.5, 0.00589, -2460600.0, -589.2] ,
            'a_0': 1.8e-05 ,
            'K_0': 2.03e+11 ,
            'Kprime_0': 4.0 ,
            'Kdprime_0': -2e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class cpv (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca1.0Si1.0O3.0'
       self.params = {
            'name': 'cpv',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1541730.0 ,
            'S_0': 73.5 ,
            'V_0': 2.745e-05 ,
            'Cp': [159.3, 0.0, -967300.0, -1075.4] ,
            'a_0': 1.87e-05 ,
            'K_0': 2.36e+11 ,
            'Kprime_0': 3.9 ,
            'Kdprime_0': -1.6e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class mak (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg1.0Si1.0O3.0'
       self.params = {
            'name': 'mak',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1490870.0 ,
            'S_0': 59.3 ,
            'V_0': 2.635e-05 ,
            'Cp': [147.8, 0.002015, -2395000.0, -801.8] ,
            'a_0': 2.12e-05 ,
            'K_0': 2.11e+11 ,
            'Kprime_0': 4.55 ,
            'Kdprime_0': -2.2e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class fak (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe1.0Si1.0O3.0'
       self.params = {
            'name': 'fak',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1142140.0 ,
            'S_0': 91.5 ,
            'V_0': 2.76e-05 ,
            'Cp': [100.3, 0.013328, -4364900.0, 419.8] ,
            'a_0': 2.12e-05 ,
            'K_0': 2.18e+11 ,
            'Kprime_0': 4.55 ,
            'Kdprime_0': -2.2e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class maj (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg4.0Si4.0O12.0'
       self.params = {
            'name': 'maj',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -6050400.0 ,
            'S_0': 255.2 ,
            'V_0': 0.00011457 ,
            'Cp': [713.6, -0.000997, -1158200.0, -6622.3] ,
            'a_0': 1.83e-05 ,
            'K_0': 1.6e+11 ,
            'Kprime_0': 4.56 ,
            'Kdprime_0': -2.8e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class py (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg3.0Al2.0Si3.0O12.0'
       self.params = {
            'name': 'py',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -6281960.0 ,
            'S_0': 269.5 ,
            'V_0': 0.00011313 ,
            'Cp': [633.5, 0.0, -5196100.0, -4315.2] ,
            'a_0': 2.37e-05 ,
            'K_0': 1.743e+11 ,
            'Kprime_0': 4.05 ,
            'Kdprime_0': -2.3e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class alm (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe3.0Al2.0Si3.0O12.0'
       self.params = {
            'name': 'alm',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -5260700.0 ,
            'S_0': 342.0 ,
            'V_0': 0.00011525 ,
            'Cp': [677.3, 0.0, -3772700.0, -5044.0] ,
            'a_0': 2.12e-05 ,
            'K_0': 1.9e+11 ,
            'Kprime_0': 2.98 ,
            'Kdprime_0': -1.6e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class spss (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mn3.0Al2.0Si3.0O12.0'
       self.params = {
            'name': 'spss',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -5693490.0 ,
            'S_0': 335.3 ,
            'V_0': 0.00011792 ,
            'Cp': [646.9, 0.0, -4525800.0, -4452.8] ,
            'a_0': 2.27e-05 ,
            'K_0': 1.74e+11 ,
            'Kprime_0': 6.68 ,
            'Kdprime_0': -3.8e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class gr (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca3.0Al2.0Si3.0O12.0'
       self.params = {
            'name': 'gr',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -6643010.0 ,
            'S_0': 255.0 ,
            'V_0': 0.00012535 ,
            'Cp': [626.0, 0.0, -5779200.0, -4002.9] ,
            'a_0': 2.2e-05 ,
            'K_0': 1.72e+11 ,
            'Kprime_0': 5.53 ,
            'Kdprime_0': -3.2e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class andr (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca3.0Fe2.0Si3.0O12.0'
       self.params = {
            'name': 'andr',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -5769100.0 ,
            'S_0': 316.4 ,
            'V_0': 0.00013204 ,
            'Cp': [638.6, 0.0, -4955100.0, -3989.2] ,
            'a_0': 2.86e-05 ,
            'K_0': 1.588e+11 ,
            'Kprime_0': 5.68 ,
            'Kdprime_0': -3.6e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class knor (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg3.0Cr2.0Si3.0O12.0'
       self.params = {
            'name': 'knor',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -5687710.0 ,
            'S_0': 317.0 ,
            'V_0': 0.00011738 ,
            'Cp': [613.0, 0.003606, -4178000.0, -3729.4] ,
            'a_0': 2.37e-05 ,
            'K_0': 1.743e+11 ,
            'Kprime_0': 4.05 ,
            'Kdprime_0': -2.3e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class osma (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='K1.0Mg2.0Al5.0Si10.0O30.0'
       self.params = {
            'name': 'osma',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -14896310.0 ,
            'S_0': 755.0 ,
            'V_0': 0.00037893 ,
            'Cp': [1540.7, -0.011359, -10339000.0, -11699.0] ,
            'a_0': 4.7e-06 ,
            'K_0': 1.29e+11 ,
            'Kprime_0': 4.1 ,
            'Kdprime_0': -3.1e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class osmm (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='K1.0Mg3.0Al3.0Si11.0O30.0'
       self.params = {
            'name': 'osmm',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -14786740.0 ,
            'S_0': 740.0 ,
            'V_0': 0.0003844 ,
            'Cp': [1525.5, -0.010267, -10538000.0, -11337.0] ,
            'a_0': 4.7e-06 ,
            'K_0': 1.29e+11 ,
            'Kprime_0': 4.1 ,
            'Kdprime_0': -3.1e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class osfa (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='K1.0Fe2.0Al5.0Si10.0O30.0'
       self.params = {
            'name': 'osfa',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -14215490.0 ,
            'S_0': 780.0 ,
            'V_0': 0.0003845 ,
            'Cp': [1558.6, -0.011359, -9476500.0, -11845.0] ,
            'a_0': 4.9e-06 ,
            'K_0': 1.29e+11 ,
            'Kprime_0': 4.1 ,
            'Kdprime_0': -3.1e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class vsv (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca19.0Mg2.0Al11.0Si18.0O78.0H9.0'
       self.params = {
            'name': 'vsv',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -42345820.0 ,
            'S_0': 1890.0 ,
            'V_0': 0.000852 ,
            'Cp': [4488.0, -0.057952, -22269300.0, -33478.0] ,
            'a_0': 2.75e-05 ,
            'K_0': 1.255e+11 ,
            'Kprime_0': 4.8 ,
            'Kdprime_0': -3.8e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class andalusite (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Al2.0Si1.0O5.0'
       self.params = {
            'name': 'and',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2588670.0 ,
            'S_0': 92.7 ,
            'V_0': 5.153e-05 ,
            'Cp': [277.3, -0.006588, -1914100.0, -2265.6] ,
            'a_0': 1.81e-05 ,
            'K_0': 1.442e+11 ,
            'Kprime_0': 6.89 ,
            'Kdprime_0': -4.8e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class ky (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Al2.0Si1.0O5.0'
       self.params = {
            'name': 'ky',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2592970.0 ,
            'S_0': 83.5 ,
            'V_0': 4.414e-05 ,
            'Cp': [279.4, -0.007124, -2055600.0, -2289.4] ,
            'a_0': 1.92e-05 ,
            'K_0': 1.601e+11 ,
            'Kprime_0': 4.05 ,
            'Kdprime_0': -2.5e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class sill (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Al2.0Si1.0O5.0'
       self.params = {
            'name': 'sill',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2585790.0 ,
            'S_0': 95.4 ,
            'V_0': 4.986e-05 ,
            'Cp': [280.2, -0.0069, -1375700.0, -2399.4] ,
            'a_0': 1.12e-05 ,
            'K_0': 1.64e+11 ,
            'Kprime_0': 5.06 ,
            'Kdprime_0': -3.1e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'BW_deltaH': 4750.0 ,
            'BW_deltaV': 1e-07 ,
            'BW_W': 4750.0 ,
            'BW_Wv': 1e-07 ,
            'BW_n': 1.0 ,
            'BW_factor': 0.25 }

class smul (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Al2.0Si1.0O5.0'
       self.params = {
            'name': 'smul',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2569210.0 ,
            'S_0': 101.5 ,
            'V_0': 4.987e-05 ,
            'Cp': [280.2, -0.0069, -1375700.0, -2399.4] ,
            'a_0': 1.36e-05 ,
            'K_0': 1.74e+11 ,
            'Kprime_0': 4.0 ,
            'Kdprime_0': -2.3e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class amul (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Al2.5Si0.5O4.75'
       self.params = {
            'name': 'amul',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2485530.0 ,
            'S_0': 113.0 ,
            'V_0': 5.083e-05 ,
            'Cp': [244.8, 0.000968, -2533300.0, -1641.6] ,
            'a_0': 1.36e-05 ,
            'K_0': 1.74e+11 ,
            'Kprime_0': 4.0 ,
            'Kdprime_0': -2.3e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class tpz (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Al2.0Si1.0O6.0H2.0'
       self.params = {
            'name': 'tpz',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2900710.0 ,
            'S_0': 100.5 ,
            'V_0': 5.339e-05 ,
            'Cp': [387.7, -0.00712, -857200.0, -3744.2] ,
            'a_0': 1.57e-05 ,
            'K_0': 1.315e+11 ,
            'Kprime_0': 4.06 ,
            'Kdprime_0': -3.1e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class mst (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg4.0Al18.0Si7.5O48.0H4.0'
       self.params = {
            'name': 'mst',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -25123740.0 ,
            'S_0': 910.0 ,
            'V_0': 0.0004426 ,
            'Cp': [2820.5, -0.059366, -13774000.0, -24126.0] ,
            'a_0': 1.81e-05 ,
            'K_0': 1.684e+11 ,
            'Kprime_0': 4.05 ,
            'Kdprime_0': -2.4e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class fst (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe4.0Al18.0Si7.5O48.0H4.0'
       self.params = {
            'name': 'fst',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -23754630.0 ,
            'S_0': 1010.0 ,
            'V_0': 0.0004488 ,
            'Cp': [2880.0, -0.056595, -10642000.0, -25373.0] ,
            'a_0': 1.83e-05 ,
            'K_0': 1.8e+11 ,
            'Kprime_0': 4.76 ,
            'Kdprime_0': -2.6e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class mnst (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mn4.0Al18.0Si7.5O48.0H4.0'
       self.params = {
            'name': 'mnst',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -24245850.0 ,
            'S_0': 1034.0 ,
            'V_0': 0.0004546 ,
            'Cp': [2873.3, -0.089064, -12688000.0, -24749.0] ,
            'a_0': 2.09e-05 ,
            'K_0': 1.8e+11 ,
            'Kprime_0': 4.76 ,
            'Kdprime_0': -2.6e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class mctd (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg1.0Al2.0Si1.0O7.0H2.0'
       self.params = {
            'name': 'mctd',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -3549250.0 ,
            'S_0': 146.0 ,
            'V_0': 6.875e-05 ,
            'Cp': [417.4, -0.003771, -2920600.0, -3417.8] ,
            'a_0': 2.63e-05 ,
            'K_0': 1.456e+11 ,
            'Kprime_0': 4.06 ,
            'Kdprime_0': -2.8e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class fctd (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe1.0Al2.0Si1.0O7.0H2.0'
       self.params = {
            'name': 'fctd',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -3208290.0 ,
            'S_0': 167.0 ,
            'V_0': 6.98e-05 ,
            'Cp': [416.1, -0.003477, -2835900.0, -3360.3] ,
            'a_0': 2.8e-05 ,
            'K_0': 1.456e+11 ,
            'Kprime_0': 4.06 ,
            'Kdprime_0': -2.8e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class mnctd (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mn1.0Al2.0Si1.0O7.0H2.0'
       self.params = {
            'name': 'mnctd',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -3336150.0 ,
            'S_0': 166.0 ,
            'V_0': 7.175e-05 ,
            'Cp': [464.4, -0.012654, -1147200.0, -4341.0] ,
            'a_0': 2.6e-05 ,
            'K_0': 1.456e+11 ,
            'Kprime_0': 4.06 ,
            'Kdprime_0': -2.8e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class merw (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca3.0Mg1.0Si2.0O8.0'
       self.params = {
            'name': 'merw',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -4545700.0 ,
            'S_0': 253.1 ,
            'V_0': 9.847e-05 ,
            'Cp': [417.5, 0.008117, -2923000.0, -2320.3] ,
            'a_0': 3.19e-05 ,
            'K_0': 1.2e+11 ,
            'Kprime_0': 4.07 ,
            'Kdprime_0': -3.4e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class spu (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca5.0Si2.0C1.0O11.0'
       self.params = {
            'name': 'spu',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -5846720.0 ,
            'S_0': 332.0 ,
            'V_0': 0.00014697 ,
            'Cp': [614.1, -0.003508, -2493100.0, -4168.0] ,
            'a_0': 3.4e-05 ,
            'K_0': 95000000000.0 ,
            'Kprime_0': 4.09 ,
            'Kdprime_0': -4.3e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class zo (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca2.0Al3.0Si3.0O13.0H1.0'
       self.params = {
            'name': 'zo',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -6896290.0 ,
            'S_0': 298.0 ,
            'V_0': 0.00013575 ,
            'Cp': [662.0, 0.010416, -6006400.0, -4260.7] ,
            'a_0': 3.12e-05 ,
            'K_0': 1.044e+11 ,
            'Kprime_0': 4.0 ,
            'Kdprime_0': -3.8e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class cz (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca2.0Al3.0Si3.0O13.0H1.0'
       self.params = {
            'name': 'cz',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -6895540.0 ,
            'S_0': 301.0 ,
            'V_0': 0.0001363 ,
            'Cp': [630.9, 0.013693, -6645800.0, -3731.1] ,
            'a_0': 2.33e-05 ,
            'K_0': 1.197e+11 ,
            'Kprime_0': 4.07 ,
            'Kdprime_0': -3.4e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class ep (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca2.0Al2.0Fe1.0Si3.0O13.0H1.0'
       self.params = {
            'name': 'ep',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -6473830.0 ,
            'S_0': 315.0 ,
            'V_0': 0.0001392 ,
            'Cp': [613.3, 0.02207, -7160000.0, -2987.7] ,
            'a_0': 2.34e-05 ,
            'K_0': 1.34e+11 ,
            'Kprime_0': 4.0 ,
            'Kdprime_0': -3e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class fep (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca2.0Al1.0Fe2.0Si3.0O13.0H1.0'
       self.params = {
            'name': 'fep',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -6028590.0 ,
            'S_0': 329.0 ,
            'V_0': 0.0001421 ,
            'Cp': [584.7, 0.030447, -7674200.0, -2244.3] ,
            'a_0': 2.31e-05 ,
            'K_0': 1.513e+11 ,
            'Kprime_0': 4.0 ,
            'Kdprime_0': -2.6e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class pmt (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca2.0Al2.0Mn1.0Si3.0O13.0H1.0'
       self.params = {
            'name': 'pmt',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -6543030.0 ,
            'S_0': 340.0 ,
            'V_0': 0.0001382 ,
            'Cp': [569.8, 0.02779, -5442900.0, -2812.6] ,
            'a_0': 2.38e-05 ,
            'K_0': 1.197e+11 ,
            'Kprime_0': 4.07 ,
            'Kdprime_0': -3.4e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class law (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca1.0Al2.0Si2.0O10.0H4.0'
       self.params = {
            'name': 'law',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -4868630.0 ,
            'S_0': 229.0 ,
            'V_0': 0.00010132 ,
            'Cp': [687.8, 0.001566, 375900.0, -7179.2] ,
            'a_0': 2.65e-05 ,
            'K_0': 1.229e+11 ,
            'Kprime_0': 5.45 ,
            'Kdprime_0': -4.4e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class mpm (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca4.0Al5.0Mg1.0Si6.0O28.0H7.0'
       self.params = {
            'name': 'mpm',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -14386910.0 ,
            'S_0': 629.0 ,
            'V_0': 0.0002955 ,
            'Cp': [1720.8, -0.024928, -5998700.0, -14620.3] ,
            'a_0': 2.48e-05 ,
            'K_0': 1.615e+11 ,
            'Kprime_0': 4.05 ,
            'Kdprime_0': -2.5e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class fpm (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca4.0Al5.0Fe1.0Si6.0O28.0H7.0'
       self.params = {
            'name': 'fpm',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -14034040.0 ,
            'S_0': 657.0 ,
            'V_0': 0.0002968 ,
            'Cp': [1737.2, -0.024582, -5161100.0, -14963.0] ,
            'a_0': 2.49e-05 ,
            'K_0': 1.615e+11 ,
            'Kprime_0': 4.05 ,
            'Kdprime_0': -2.5e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class jgd (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca4.0Fe6.0Si6.0O28.0H7.0'
       self.params = {
            'name': 'jgd',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -11808960.0 ,
            'S_0': 830.0 ,
            'V_0': 0.0003108 ,
            'Cp': [1795.4, -0.037986, -4455700.0, -14888.0] ,
            'a_0': 2.49e-05 ,
            'K_0': 1.615e+11 ,
            'Kprime_0': 4.05 ,
            'Kdprime_0': -2.5e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class geh (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca2.0Al2.0Si1.0O7.0'
       self.params = {
            'name': 'geh',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -3992240.0 ,
            'S_0': 198.5 ,
            'V_0': 9.024e-05 ,
            'Cp': [405.7, -0.007099, -1188300.0, -3174.4] ,
            'a_0': 2.23e-05 ,
            'K_0': 1.08e+11 ,
            'Kprime_0': 4.08 ,
            'Kdprime_0': -3.8e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'BW_deltaH': 7510.0 ,
            'BW_deltaV': 9e-07 ,
            'BW_W': 7500.0 ,
            'BW_Wv': 9e-07 ,
            'BW_n': 1.0 ,
            'BW_factor': 0.8 }

class ak (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca2.0Mg1.0Si2.0O7.0'
       self.params = {
            'name': 'ak',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -3865620.0 ,
            'S_0': 212.5 ,
            'V_0': 9.254e-05 ,
            'Cp': [385.4, 0.003209, -247500.0, -2889.9] ,
            'a_0': 2.57e-05 ,
            'K_0': 1.42e+11 ,
            'Kprime_0': 4.06 ,
            'Kdprime_0': -2.9e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class rnk (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca3.0Si2.0O7.0'
       self.params = {
            'name': 'rnk',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -3943820.0 ,
            'S_0': 210.0 ,
            'V_0': 9.651e-05 ,
            'Cp': [372.3, -0.002893, -2462400.0, -2181.3] ,
            'a_0': 3.28e-05 ,
            'K_0': 95000000000.0 ,
            'Kprime_0': 4.09 ,
            'Kdprime_0': -4.3e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class ty (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca5.0Si2.0C2.0O13.0'
       self.params = {
            'name': 'ty',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -6368040.0 ,
            'S_0': 390.0 ,
            'V_0': 0.00017039 ,
            'Cp': [741.7, -0.005345, -1434600.0, -5878.5] ,
            'a_0': 3.42e-05 ,
            'K_0': 95000000000.0 ,
            'Kprime_0': 4.09 ,
            'Kdprime_0': -4.3e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class crd (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg2.0Al4.0Si5.0O18.0'
       self.params = {
            'name': 'crd',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -9163430.0 ,
            'S_0': 404.1 ,
            'V_0': 0.00023322 ,
            'Cp': [906.1, 0.0, -7902000.0, -6293.4] ,
            'a_0': 6.8e-06 ,
            'K_0': 1.29e+11 ,
            'Kprime_0': 4.1 ,
            'Kdprime_0': -3.1e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'BW_deltaH': 36710.0 ,
            'BW_deltaV': 1e-06 ,
            'BW_W': 36700.0 ,
            'BW_Wv': 1e-06 ,
            'BW_n': 2.0 ,
            'BW_factor': 1.5 }

class hcrd (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg2.0Al4.0Si5.0O19.0H2.0'
       self.params = {
            'name': 'hcrd',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -9448520.0 ,
            'S_0': 483.0 ,
            'V_0': 0.00023322 ,
            'Cp': [955.3, 0.0, -8352600.0, -6301.2] ,
            'a_0': 6.7e-06 ,
            'K_0': 1.29e+11 ,
            'Kprime_0': 4.1 ,
            'Kdprime_0': -3.1e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'BW_deltaH': 36710.0 ,
            'BW_deltaV': 1e-06 ,
            'BW_W': 36700.0 ,
            'BW_Wv': 1e-06 ,
            'BW_n': 2.0 ,
            'BW_factor': 1.5 }

class fcrd (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe2.0Al4.0Si5.0O18.0'
       self.params = {
            'name': 'fcrd',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -8444070.0 ,
            'S_0': 461.0 ,
            'V_0': 0.0002371 ,
            'Cp': [924.0, 0.0, -7039400.0, -6439.6] ,
            'a_0': 6.7e-06 ,
            'K_0': 1.29e+11 ,
            'Kprime_0': 4.1 ,
            'Kdprime_0': -3.1e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'BW_deltaH': 36710.0 ,
            'BW_deltaV': 1e-06 ,
            'BW_W': 36700.0 ,
            'BW_Wv': 1e-06 ,
            'BW_n': 2.0 ,
            'BW_factor': 1.5 }

class mncrd (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mn2.0Al4.0Si5.0O18.0'
       self.params = {
            'name': 'mncrd',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -8693590.0 ,
            'S_0': 473.0 ,
            'V_0': 0.00024027 ,
            'Cp': [886.5, 0.0, -8840000.0, -5590.4] ,
            'a_0': 6.9e-06 ,
            'K_0': 1.29e+11 ,
            'Kprime_0': 4.1 ,
            'Kdprime_0': -3.1e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'BW_deltaH': 36710.0 ,
            'BW_deltaV': 1e-06 ,
            'BW_W': 36700.0 ,
            'BW_Wv': 1e-06 ,
            'BW_n': 2.0 ,
            'BW_factor': 1.5 }

class phA (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg7.0Si2.0O14.0H6.0'
       self.params = {
            'name': 'phA',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -7129620.0 ,
            'S_0': 350.5 ,
            'V_0': 0.00015422 ,
            'Cp': [962.0, -0.011521, -4517800.0, -7724.7] ,
            'a_0': 3.55e-05 ,
            'K_0': 1.45e+11 ,
            'Kprime_0': 4.06 ,
            'Kdprime_0': -2.8e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class sph (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca1.0Ti1.0Si1.0O5.0'
       self.params = {
            'name': 'sph',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2601660.0 ,
            'S_0': 124.0 ,
            'V_0': 5.565e-05 ,
            'Cp': [227.9, 0.002924, -3539500.0, -894.3] ,
            'a_0': 1.58e-05 ,
            'K_0': 1.017e+11 ,
            'Kprime_0': 9.85 ,
            'Kdprime_0': -9.7e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'landau_Tc': 485.0 ,
            'landau_Smax': 0.4 ,
            'landau_Vmax': 5e-08 }

class cstn (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca1.0Si2.0O5.0'
       self.params = {
            'name': 'cstn',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2496350.0 ,
            'S_0': 99.5 ,
            'V_0': 4.818e-05 ,
            'Cp': [205.6, 0.006034, -5517700.0, -352.6] ,
            'a_0': 1.58e-05 ,
            'K_0': 1.782e+11 ,
            'Kprime_0': 4.0 ,
            'Kdprime_0': -2.2e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class zrc (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Zr1.0Si1.0O4.0'
       self.params = {
            'name': 'zrc',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2035070.0 ,
            'S_0': 83.03 ,
            'V_0': 3.926e-05 ,
            'Cp': [232.0, -0.014405, 0.0, -2238.2] ,
            'a_0': 1.25e-05 ,
            'K_0': 2.301e+11 ,
            'Kprime_0': 4.04 ,
            'Kdprime_0': -1.8e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class en (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg2.0Si2.0O6.0'
       self.params = {
            'name': 'en',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -3090220.0 ,
            'S_0': 132.5 ,
            'V_0': 6.262e-05 ,
            'Cp': [356.2, -0.00299, -596900.0, -3185.3] ,
            'a_0': 2.27e-05 ,
            'K_0': 1.059e+11 ,
            'Kprime_0': 8.65 ,
            'Kdprime_0': -8.2e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class pren (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg2.0Si2.0O6.0'
       self.params = {
            'name': 'pren',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -3084560.0 ,
            'S_0': 137.0 ,
            'V_0': 6.476e-05 ,
            'Cp': [356.2, -0.00299, -596900.0, -3185.3] ,
            'a_0': 2.3e-05 ,
            'K_0': 1.059e+11 ,
            'Kprime_0': 8.65 ,
            'Kdprime_0': -8.2e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class cen (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg2.0Si2.0O6.0'
       self.params = {
            'name': 'cen',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -3091110.0 ,
            'S_0': 132.0 ,
            'V_0': 6.264e-05 ,
            'Cp': [306.0, -0.003793, -3041700.0, -1852.1] ,
            'a_0': 2.11e-05 ,
            'K_0': 1.059e+11 ,
            'Kprime_0': 8.65 ,
            'Kdprime_0': -8.2e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class hen (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg2.0Si2.0O6.0'
       self.params = {
            'name': 'hen',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -3082730.0 ,
            'S_0': 131.7 ,
            'V_0': 6.099e-05 ,
            'Cp': [356.2, -0.00299, -596900.0, -3185.3] ,
            'a_0': 2.26e-05 ,
            'K_0': 1.5e+11 ,
            'Kprime_0': 5.5 ,
            'Kdprime_0': -3.6e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class fs (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe2.0Si2.0O6.0'
       self.params = {
            'name': 'fs',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2388710.0 ,
            'S_0': 189.9 ,
            'V_0': 6.592e-05 ,
            'Cp': [398.7, -0.006579, 1290100.0, -4058.0] ,
            'a_0': 3.26e-05 ,
            'K_0': 1.01e+11 ,
            'Kprime_0': 4.08 ,
            'Kdprime_0': -4e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class mgts (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg1.0Al2.0Si1.0O6.0'
       self.params = {
            'name': 'mgts',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -3196670.0 ,
            'S_0': 131.0 ,
            'V_0': 6.05e-05 ,
            'Cp': [371.4, -0.004082, -398400.0, -3547.1] ,
            'a_0': 2.17e-05 ,
            'K_0': 1.028e+11 ,
            'Kprime_0': 8.55 ,
            'Kdprime_0': -8.3e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class di (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca1.0Mg1.0Si2.0O6.0'
       self.params = {
            'name': 'di',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -3201850.0 ,
            'S_0': 142.9 ,
            'V_0': 6.619e-05 ,
            'Cp': [314.5, 4.1e-05, -2745900.0, -2020.1] ,
            'a_0': 2.73e-05 ,
            'K_0': 1.192e+11 ,
            'Kprime_0': 5.19 ,
            'Kdprime_0': -4.4e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class hed (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca1.0Fe1.0Si2.0O6.0'
       self.params = {
            'name': 'hed',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2842060.0 ,
            'S_0': 175.0 ,
            'V_0': 6.795e-05 ,
            'Cp': [340.2, 0.000812, -1047800.0, -2646.7] ,
            'a_0': 2.38e-05 ,
            'K_0': 1.192e+11 ,
            'Kprime_0': 3.97 ,
            'Kdprime_0': -3.3e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class jd (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Na1.0Al1.0Si2.0O6.0'
       self.params = {
            'name': 'jd',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -3025270.0 ,
            'S_0': 133.5 ,
            'V_0': 6.04e-05 ,
            'Cp': [319.4, 0.003616, -1173900.0, -2469.5] ,
            'a_0': 2.1e-05 ,
            'K_0': 1.281e+11 ,
            'Kprime_0': 3.81 ,
            'Kdprime_0': -3e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class acm (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Na1.0Fe1.0Si2.0O6.0'
       self.params = {
            'name': 'acm',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2583430.0 ,
            'S_0': 170.6 ,
            'V_0': 6.459e-05 ,
            'Cp': [307.1, 0.016758, -1685500.0, -2125.8] ,
            'a_0': 2.11e-05 ,
            'K_0': 1.06e+11 ,
            'Kprime_0': 4.08 ,
            'Kdprime_0': -3.8e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class kos (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Na1.0Cr1.0Si2.0O6.0'
       self.params = {
            'name': 'kos',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2746840.0 ,
            'S_0': 149.65 ,
            'V_0': 6.309e-05 ,
            'Cp': [309.2, 0.005419, -664600.0, -2176.6] ,
            'a_0': 1.94e-05 ,
            'K_0': 1.308e+11 ,
            'Kprime_0': 3.0 ,
            'Kdprime_0': -2.3e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class cats (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca1.0Al2.0Si1.0O6.0'
       self.params = {
            'name': 'cats',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -3310110.0 ,
            'S_0': 135.0 ,
            'V_0': 6.356e-05 ,
            'Cp': [347.6, -0.006974, -1781600.0, -2757.5] ,
            'a_0': 2.08e-05 ,
            'K_0': 1.192e+11 ,
            'Kprime_0': 5.19 ,
            'Kdprime_0': -4.4e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'BW_deltaH': 3800.0 ,
            'BW_deltaV': 1e-07 ,
            'BW_W': 3800.0 ,
            'BW_Wv': 1e-07 ,
            'BW_n': 1.0 ,
            'BW_factor': 0.25 }

class caes (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca0.5Al1.0Si2.0O6.0'
       self.params = {
            'name': 'caes',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -3002020.0 ,
            'S_0': 127.0 ,
            'V_0': 6.05e-05 ,
            'Cp': [362.0, -0.016944, -175900.0, -3565.7] ,
            'a_0': 2.31e-05 ,
            'K_0': 1.192e+11 ,
            'Kprime_0': 5.19 ,
            'Kdprime_0': -4.4e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class rhod (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mn1.0Si1.0O3.0'
       self.params = {
            'name': 'rhod',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1322380.0 ,
            'S_0': 100.5 ,
            'V_0': 3.494e-05 ,
            'Cp': [138.4, 0.004088, -1936000.0, -538.9] ,
            'a_0': 2.81e-05 ,
            'K_0': 84000000000.0 ,
            'Kprime_0': 4.0 ,
            'Kdprime_0': -4.8e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class pxmn (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mn1.0Si1.0O3.0'
       self.params = {
            'name': 'pxmn',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1323160.0 ,
            'S_0': 99.3 ,
            'V_0': 3.472e-05 ,
            'Cp': [138.4, 0.004088, -1936000.0, -538.9] ,
            'a_0': 2.8e-05 ,
            'K_0': 84000000000.0 ,
            'Kprime_0': 4.0 ,
            'Kdprime_0': -4.8e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class wo (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca1.0Si1.0O3.0'
       self.params = {
            'name': 'wo',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1633770.0 ,
            'S_0': 82.5 ,
            'V_0': 3.993e-05 ,
            'Cp': [159.3, 0.0, -967300.0, -1075.4] ,
            'a_0': 2.54e-05 ,
            'K_0': 79500000000.0 ,
            'Kprime_0': 4.1 ,
            'Kdprime_0': -5.2e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class pswo (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca1.0Si1.0O3.0'
       self.params = {
            'name': 'pswo',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1627960.0 ,
            'S_0': 87.8 ,
            'V_0': 4.008e-05 ,
            'Cp': [157.8, 0.0, -967300.0, -1075.4] ,
            'a_0': 2.85e-05 ,
            'K_0': 1.1e+11 ,
            'Kprime_0': 4.08 ,
            'Kdprime_0': -3.7e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class wal (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca1.0Si1.0O3.0'
       self.params = {
            'name': 'wal',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1625900.0 ,
            'S_0': 83.5 ,
            'V_0': 3.7633e-05 ,
            'Cp': [159.3, 0.0, -967300.0, -1075.4] ,
            'a_0': 2.54e-05 ,
            'K_0': 79500000000.0 ,
            'Kprime_0': 4.1 ,
            'Kdprime_0': -5.2e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class tr (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca2.0Mg5.0Si8.0O24.0H2.0'
       self.params = {
            'name': 'tr',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -12304870.0 ,
            'S_0': 553.0 ,
            'V_0': 0.0002727 ,
            'Cp': [1260.2, 0.00383, -11455000.0, -8237.6] ,
            'a_0': 2.61e-05 ,
            'K_0': 76200000000.0 ,
            'Kprime_0': 4.1 ,
            'Kdprime_0': -5.4e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class fact (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca2.0Fe5.0Si8.0O24.0H2.0'
       self.params = {
            'name': 'fact',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -10504120.0 ,
            'S_0': 710.0 ,
            'V_0': 0.0002842 ,
            'Cp': [1290.0, 0.029992, -8447500.0, -8947.0] ,
            'a_0': 2.88e-05 ,
            'K_0': 76000000000.0 ,
            'Kprime_0': 4.1 ,
            'Kdprime_0': -5.4e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class ts (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca2.0Mg3.0Al4.0Si6.0O24.0H2.0'
       self.params = {
            'name': 'ts',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -12555270.0 ,
            'S_0': 533.0 ,
            'V_0': 0.000268 ,
            'Cp': [1244.8, 0.024348, -11965000.0, -8112.1] ,
            'a_0': 2.66e-05 ,
            'K_0': 76000000000.0 ,
            'Kprime_0': 4.1 ,
            'Kdprime_0': -5.4e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class parg (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Na1.0Ca2.0Mg4.0Al3.0Si6.0O24.0H2.0'
       self.params = {
            'name': 'parg',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -12664730.0 ,
            'S_0': 635.0 ,
            'V_0': 0.0002719 ,
            'Cp': [1280.2, 0.022997, -12359500.0, -8065.8] ,
            'a_0': 2.8e-05 ,
            'K_0': 91200000000.0 ,
            'Kprime_0': 4.09 ,
            'Kdprime_0': -4.5e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class gl (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Na2.0Mg3.0Al2.0Si8.0O24.0H2.0'
       self.params = {
            'name': 'gl',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -11960240.0 ,
            'S_0': 530.0 ,
            'V_0': 0.0002598 ,
            'Cp': [1717.5, -0.12107, 7075000.0, -19272.0] ,
            'a_0': 1.49e-05 ,
            'K_0': 88300000000.0 ,
            'Kprime_0': 4.09 ,
            'Kdprime_0': -4.6e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class fgl (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Na2.0Al2.0Fe3.0Si8.0O24.0H2.0'
       self.params = {
            'name': 'fgl',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -10880210.0 ,
            'S_0': 624.0 ,
            'V_0': 0.0002659 ,
            'Cp': [1762.9, -0.118992, 9423700.0, -20207.1] ,
            'a_0': 1.83e-05 ,
            'K_0': 89000000000.0 ,
            'Kprime_0': 4.09 ,
            'Kdprime_0': -4.6e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class rieb (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Na2.0Fe5.0Si8.0O24.0H2.0'
       self.params = {
            'name': 'rieb',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -10024780.0 ,
            'S_0': 695.0 ,
            'V_0': 0.0002749 ,
            'Cp': [1787.3, -0.124882, 9627100.0, -20275.5] ,
            'a_0': 1.81e-05 ,
            'K_0': 89000000000.0 ,
            'Kprime_0': 4.09 ,
            'Kdprime_0': -4.6e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class anth (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg7.0Si8.0O24.0H2.0'
       self.params = {
            'name': 'anth',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -12066840.0 ,
            'S_0': 537.0 ,
            'V_0': 0.0002654 ,
            'Cp': [1277.3, 0.025825, -9704600.0, -9074.7] ,
            'a_0': 2.52e-05 ,
            'K_0': 70000000000.0 ,
            'Kprime_0': 4.11 ,
            'Kdprime_0': -5.9e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class fanth (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe7.0Si8.0O24.0H2.0'
       self.params = {
            'name': 'fanth',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -9624520.0 ,
            'S_0': 725.0 ,
            'V_0': 0.0002787 ,
            'Cp': [1383.1, 0.030669, -4224700.0, -11257.6] ,
            'a_0': 2.74e-05 ,
            'K_0': 70000000000.0 ,
            'Kprime_0': 4.11 ,
            'Kdprime_0': -5.9e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class cumm (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg7.0Si8.0O24.0H2.0'
       self.params = {
            'name': 'cumm',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -12064690.0 ,
            'S_0': 538.0 ,
            'V_0': 0.0002633 ,
            'Cp': [1277.3, 0.025825, -9704600.0, -9074.7] ,
            'a_0': 2.52e-05 ,
            'K_0': 70000000000.0 ,
            'Kprime_0': 4.11 ,
            'Kdprime_0': -5.9e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class grun (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe7.0Si8.0O24.0H2.0'
       self.params = {
            'name': 'grun',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -9607150.0 ,
            'S_0': 735.0 ,
            'V_0': 0.0002784 ,
            'Cp': [1383.1, 0.030669, -4224700.0, -11257.6] ,
            'a_0': 2.74e-05 ,
            'K_0': 64800000000.0 ,
            'Kprime_0': 4.12 ,
            'Kdprime_0': -6.4e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class ged (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg5.0Al4.0Si6.0O24.0H2.0'
       self.params = {
            'name': 'ged',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -12329140.0 ,
            'S_0': 517.0 ,
            'V_0': 0.00025548 ,
            'Cp': [1307.7, 0.023642, -9307400.0, -9799.0] ,
            'a_0': 2.41e-05 ,
            'K_0': 77000000000.0 ,
            'Kprime_0': 4.1 ,
            'Kdprime_0': -5.3e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class spr4 (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg4.0Al8.0Si2.0O20.0'
       self.params = {
            'name': 'spr4',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -11022020.0 ,
            'S_0': 425.5 ,
            'V_0': 0.000199 ,
            'Cp': [1133.1, -0.007596, -8816600.0, -8180.6] ,
            'a_0': 2.05e-05 ,
            'K_0': 2.5e+11 ,
            'Kprime_0': 4.04 ,
            'Kdprime_0': -1.6e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class spr5 (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg3.0Al10.0Si1.0O20.0'
       self.params = {
            'name': 'spr5',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -11135570.0 ,
            'S_0': 419.5 ,
            'V_0': 0.0001975 ,
            'Cp': [1103.4, 0.001015, -10957000.0, -7409.2] ,
            'a_0': 2.06e-05 ,
            'K_0': 2.5e+11 ,
            'Kprime_0': 4.04 ,
            'Kdprime_0': -1.6e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class fspr (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe4.0Al8.0Si2.0O20.0'
       self.params = {
            'name': 'fspr',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -9659530.0 ,
            'S_0': 485.0 ,
            'V_0': 0.00019923 ,
            'Cp': [1132.9, -0.007348, -10420200.0, -7036.6] ,
            'a_0': 1.96e-05 ,
            'K_0': 2.5e+11 ,
            'Kprime_0': 4.04 ,
            'Kdprime_0': -1.7e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class mcar (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg1.0Al2.0Si2.0O10.0H4.0'
       self.params = {
            'name': 'mcar',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -4771050.0 ,
            'S_0': 221.5 ,
            'V_0': 0.0001059 ,
            'Cp': [683.0, -0.014054, 291000.0, -6976.4] ,
            'a_0': 2.43e-05 ,
            'K_0': 52500000000.0 ,
            'Kprime_0': 4.14 ,
            'Kdprime_0': -7.9e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class fcar (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe1.0Al2.0Si2.0O10.0H4.0'
       self.params = {
            'name': 'fcar',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -4411440.0 ,
            'S_0': 251.1 ,
            'V_0': 0.00010695 ,
            'Cp': [686.6, -0.012415, 186000.0, -6884.0] ,
            'a_0': 2.21e-05 ,
            'K_0': 52500000000.0 ,
            'Kprime_0': 4.14 ,
            'Kdprime_0': -7.9e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class deer (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe18.0Si12.0O50.0H10.0'
       self.params = {
            'name': 'deer',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -18341400.0 ,
            'S_0': 1650.0 ,
            'V_0': 0.0005574 ,
            'Cp': [3164.4, -0.027883, -5039100.0, -26721.0] ,
            'a_0': 2.75e-05 ,
            'K_0': 63000000000.0 ,
            'Kprime_0': 4.12 ,
            'Kdprime_0': -6.5e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class mu (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='K1.0Al3.0Si3.0O12.0H2.0'
       self.params = {
            'name': 'mu',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -5976510.0 ,
            'S_0': 292.0 ,
            'V_0': 0.00014083 ,
            'Cp': [756.4, -0.01984, -2170000.0, -6979.2] ,
            'a_0': 3.07e-05 ,
            'K_0': 49000000000.0 ,
            'Kprime_0': 4.15 ,
            'Kdprime_0': -8.5e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class cel (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='K1.0Mg1.0Al1.0Si4.0O12.0H2.0'
       self.params = {
            'name': 'cel',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -5834840.0 ,
            'S_0': 290.0 ,
            'V_0': 0.00013957 ,
            'Cp': [741.2, -0.018748, -2368800.0, -6616.9] ,
            'a_0': 3.07e-05 ,
            'K_0': 70000000000.0 ,
            'Kprime_0': 4.11 ,
            'Kdprime_0': -5.9e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class fcel (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='K1.0Fe1.0Al1.0Si4.0O12.0H2.0'
       self.params = {
            'name': 'fcel',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -5468490.0 ,
            'S_0': 330.0 ,
            'V_0': 0.0001407 ,
            'Cp': [756.3, -0.019147, -1586100.0, -6928.7] ,
            'a_0': 3.18e-05 ,
            'K_0': 70000000000.0 ,
            'Kprime_0': 4.11 ,
            'Kdprime_0': -5.9e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class pa (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Na1.0Al3.0Si3.0O12.0H2.0'
       self.params = {
            'name': 'pa',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -5942840.0 ,
            'S_0': 277.0 ,
            'V_0': 0.00013211 ,
            'Cp': [803.0, -0.03158, 217000.0, -8151.0] ,
            'a_0': 3.7e-05 ,
            'K_0': 51500000000.0 ,
            'Kprime_0': 6.51 ,
            'Kdprime_0': -1.26e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class ma (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca1.0Al4.0Si2.0O12.0H2.0'
       self.params = {
            'name': 'ma',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -6242070.0 ,
            'S_0': 265.0 ,
            'V_0': 0.00012964 ,
            'Cp': [744.4, -0.0168, -2074400.0, -6783.2] ,
            'a_0': 2.33e-05 ,
            'K_0': 1e+11 ,
            'Kprime_0': 4.08 ,
            'Kdprime_0': -4.1e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class phl (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='K1.0Mg3.0Al1.0Si3.0O12.0H2.0'
       self.params = {
            'name': 'phl',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -6214880.0 ,
            'S_0': 326.0 ,
            'V_0': 0.00014964 ,
            'Cp': [770.3, -0.036939, -2328900.0, -6531.6] ,
            'a_0': 3.8e-05 ,
            'K_0': 51300000000.0 ,
            'Kprime_0': 7.33 ,
            'Kdprime_0': -1.43e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class ann (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='K1.0Fe3.0Al1.0Si3.0O12.0H2.0'
       self.params = {
            'name': 'ann',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -5143720.0 ,
            'S_0': 420.0 ,
            'V_0': 0.00015432 ,
            'Cp': [815.7, -0.034861, 19800.0, -7466.7] ,
            'a_0': 3.8e-05 ,
            'K_0': 51300000000.0 ,
            'Kprime_0': 7.33 ,
            'Kdprime_0': -1.43e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class mnbi (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='K1.0Mn3.0Al1.0Si3.0O12.0H2.0'
       self.params = {
            'name': 'mnbi',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -5477520.0 ,
            'S_0': 433.0 ,
            'V_0': 0.00015264 ,
            'Cp': [809.9, -0.059213, -1514400.0, -6998.7] ,
            'a_0': 3.8e-05 ,
            'K_0': 53000000000.0 ,
            'Kprime_0': 7.33 ,
            'Kdprime_0': -1.43e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class east (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='K1.0Mg2.0Al3.0Si2.0O12.0H2.0'
       self.params = {
            'name': 'east',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -6330380.0 ,
            'S_0': 318.0 ,
            'V_0': 0.00014738 ,
            'Cp': [785.5, -0.038031, -2130300.0, -6893.7] ,
            'a_0': 3.8e-05 ,
            'K_0': 53000000000.0 ,
            'Kprime_0': 7.33 ,
            'Kdprime_0': -1.43e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class naph (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Na1.0Mg3.0Al1.0Si3.0O12.0H2.0'
       self.params = {
            'name': 'naph',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -6172010.0 ,
            'S_0': 318.0 ,
            'V_0': 0.0001445 ,
            'Cp': [773.5, -0.040229, -2597900.0, -6512.6] ,
            'a_0': 3.28e-05 ,
            'K_0': 51300000000.0 ,
            'Kprime_0': 7.33 ,
            'Kdprime_0': -1.43e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class clin (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg5.0Al2.0Si3.0O18.0H8.0'
       self.params = {
            'name': 'clin',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -8909160.0 ,
            'S_0': 437.0 ,
            'V_0': 0.0002114 ,
            'Cp': [1170.8, -0.001508, -3825800.0, -10315.0] ,
            'a_0': 2.04e-05 ,
            'K_0': 87000000000.0 ,
            'Kprime_0': 4.09 ,
            'Kdprime_0': -4.7e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class ames (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg4.0Al4.0Si2.0O18.0H8.0'
       self.params = {
            'name': 'ames',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -9040460.0 ,
            'S_0': 412.0 ,
            'V_0': 0.0002071 ,
            'Cp': [1186.0, -0.002599, -3627200.0, -10677.0] ,
            'a_0': 2e-05 ,
            'K_0': 87000000000.0 ,
            'Kprime_0': 4.09 ,
            'Kdprime_0': -4.7e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class afchl (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg6.0Si4.0O18.0H8.0'
       self.params = {
            'name': 'afchl',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -8727860.0 ,
            'S_0': 439.0 ,
            'V_0': 0.0002157 ,
            'Cp': [1155.0, -0.000417, -4024400.0, -9952.9] ,
            'a_0': 2.04e-05 ,
            'K_0': 87000000000.0 ,
            'Kprime_0': 4.09 ,
            'Kdprime_0': -4.7e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class daph (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe5.0Al2.0Si3.0O18.0H8.0'
       self.params = {
            'name': 'daph',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -7116910.0 ,
            'S_0': 584.0 ,
            'V_0': 0.0002162 ,
            'Cp': [1192.0, -0.00594, -4826400.0, -9768.3] ,
            'a_0': 2.27e-05 ,
            'K_0': 87000000000.0 ,
            'Kprime_0': 4.09 ,
            'Kdprime_0': -4.7e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class mnchl (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mn5.0Al2.0Si3.0O18.0H8.0'
       self.params = {
            'name': 'mnchl',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -7702320.0 ,
            'S_0': 595.0 ,
            'V_0': 0.0002259 ,
            'Cp': [1136.5, -0.005243, -5548100.0, -8911.5] ,
            'a_0': 2.23e-05 ,
            'K_0': 87000000000.0 ,
            'Kprime_0': 4.09 ,
            'Kdprime_0': -4.7e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class sud (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg2.0Al4.0Si3.0O18.0H8.0'
       self.params = {
            'name': 'sud',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -8626540.0 ,
            'S_0': 395.0 ,
            'V_0': 0.000203 ,
            'Cp': [1436.1, -0.048749, -2748500.0, -13764.0] ,
            'a_0': 1.99e-05 ,
            'K_0': 87000000000.0 ,
            'Kprime_0': 4.09 ,
            'Kdprime_0': -4.7e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class fsud (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe2.0Al4.0Si3.0O18.0H8.0'
       self.params = {
            'name': 'fsud',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -7899850.0 ,
            'S_0': 456.0 ,
            'V_0': 0.000204 ,
            'Cp': [1466.3, -0.047365, -1182800.0, -14388.0] ,
            'a_0': 2.08e-05 ,
            'K_0': 87000000000.0 ,
            'Kprime_0': 4.09 ,
            'Kdprime_0': -4.7e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class prl (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Al2.0Si4.0O12.0H2.0'
       self.params = {
            'name': 'prl',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -5640610.0 ,
            'S_0': 239.0 ,
            'V_0': 0.00012804 ,
            'Cp': [784.5, -0.042948, 1251000.0, -8495.9] ,
            'a_0': 4.5e-05 ,
            'K_0': 37000000000.0 ,
            'Kprime_0': 10.0 ,
            'Kdprime_0': -2.71e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class ta (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg3.0Si4.0O12.0H2.0'
       self.params = {
            'name': 'ta',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -5897170.0 ,
            'S_0': 259.0 ,
            'V_0': 0.00013665 ,
            'Cp': [622.2, 0.0, -6385500.0, -3916.3] ,
            'a_0': 1.8e-05 ,
            'K_0': 43000000000.0 ,
            'Kprime_0': 6.17 ,
            'Kdprime_0': -1.44e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class fta (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe3.0Si4.0O12.0H2.0'
       self.params = {
            'name': 'fta',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -4798540.0 ,
            'S_0': 352.0 ,
            'V_0': 0.00014225 ,
            'Cp': [579.7, 0.039494, -6459300.0, -3088.1] ,
            'a_0': 1.8e-05 ,
            'K_0': 43000000000.0 ,
            'Kprime_0': 6.17 ,
            'Kdprime_0': -1.44e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class tats (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg2.0Al2.0Si3.0O12.0H2.0'
       self.params = {
            'name': 'tats',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -6001290.0 ,
            'S_0': 259.0 ,
            'V_0': 0.0001351 ,
            'Cp': [549.5, 0.036324, -8606600.0, -2515.3] ,
            'a_0': 1.8e-05 ,
            'K_0': 43000000000.0 ,
            'Kprime_0': 6.17 ,
            'Kdprime_0': -1.44e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class tap (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Al2.0Si4.0O12.0H2.0'
       self.params = {
            'name': 'tap',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -5649780.0 ,
            'S_0': 235.0 ,
            'V_0': 0.0001345 ,
            'Cp': [784.5, -0.042948, 1251000.0, -8495.9] ,
            'a_0': 4.5e-05 ,
            'K_0': 37000000000.0 ,
            'Kprime_0': 10.0 ,
            'Kdprime_0': -2.71e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class minn (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe3.0Si4.0O12.0H2.0'
       self.params = {
            'name': 'minn',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -4819310.0 ,
            'S_0': 355.0 ,
            'V_0': 0.00014851 ,
            'Cp': [579.7, 0.039494, -6459300.0, -3088.1] ,
            'a_0': 1.8e-05 ,
            'K_0': 43000000000.0 ,
            'Kprime_0': 6.17 ,
            'Kdprime_0': -1.44e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class minm (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg3.0Si4.0O12.0H2.0'
       self.params = {
            'name': 'minm',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -5866000.0 ,
            'S_0': 263.9 ,
            'V_0': 0.00014291 ,
            'Cp': [622.2, 0.0, -6385500.0, -3916.3] ,
            'a_0': 1.8e-05 ,
            'K_0': 43000000000.0 ,
            'Kprime_0': 6.17 ,
            'Kdprime_0': -1.44e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class kao (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Al2.0Si2.0O9.0H4.0'
       self.params = {
            'name': 'kao',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -4122000.0 ,
            'S_0': 203.7 ,
            'V_0': 9.934e-05 ,
            'Cp': [436.7, -0.034295, -4055900.0, -2699.1] ,
            'a_0': 2.51e-05 ,
            'K_0': 64500000000.0 ,
            'Kprime_0': 4.12 ,
            'Kdprime_0': -6.4e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class pre (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca2.0Al2.0Si3.0O12.0H2.0'
       self.params = {
            'name': 'pre',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -6202170.0 ,
            'S_0': 292.8 ,
            'V_0': 0.00014026 ,
            'Cp': [724.9, -0.013865, -2059000.0, -6323.9] ,
            'a_0': 1.58e-05 ,
            'K_0': 1.093e+11 ,
            'Kprime_0': 4.01 ,
            'Kdprime_0': -3.7e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class fpre (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca2.0Al1.0Fe1.0Si3.0O12.0H2.0'
       self.params = {
            'name': 'fpre',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -5766640.0 ,
            'S_0': 320.0 ,
            'V_0': 0.000148 ,
            'Cp': [737.1, -0.01681, -1957300.0, -6358.1] ,
            'a_0': 1.58e-05 ,
            'K_0': 1.093e+11 ,
            'Kprime_0': 4.01 ,
            'Kdprime_0': -3.7e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class chr (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg3.0Si2.0O9.0H4.0'
       self.params = {
            'name': 'chr',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -4361000.0 ,
            'S_0': 221.3 ,
            'V_0': 0.00010746 ,
            'Cp': [624.7, -0.02077, -1721800.0, -5619.4] ,
            'a_0': 2.2e-05 ,
            'K_0': 62800000000.0 ,
            'Kprime_0': 4.0 ,
            'Kdprime_0': -6.4e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class liz (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg3.0Si2.0O9.0H4.0'
       self.params = {
            'name': 'liz',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -4369190.0 ,
            'S_0': 212.0 ,
            'V_0': 0.00010645 ,
            'Cp': [614.7, -0.02077, -1721800.0, -5619.4] ,
            'a_0': 2.2e-05 ,
            'K_0': 71000000000.0 ,
            'Kprime_0': 3.2 ,
            'Kdprime_0': -4.5e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class glt (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe3.0Si2.0O9.0H4.0'
       self.params = {
            'name': 'glt',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -3297620.0 ,
            'S_0': 310.0 ,
            'V_0': 0.0001198 ,
            'Cp': [576.4, 0.002984, -3757000.0, -4166.2] ,
            'a_0': 2.28e-05 ,
            'K_0': 63000000000.0 ,
            'Kprime_0': 4.0 ,
            'Kdprime_0': -6.3e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class fstp (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='K0.5Fe5.0Si8.0Al2.0O30.5H12.5'
       self.params = {
            'name': 'fstp',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -12551070.0 ,
            'S_0': 930.2 ,
            'V_0': 0.00037239 ,
            'Cp': [1944.3, -0.012289, -4840200.0, -16635.0] ,
            'a_0': 3.68e-05 ,
            'K_0': 51300000000.0 ,
            'Kprime_0': 7.33 ,
            'Kdprime_0': -1.43e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class mstp (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='K0.5Mg5.0Si8.0Al2.0O30.5H12.5'
       self.params = {
            'name': 'mstp',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -14288380.0 ,
            'S_0': 847.4 ,
            'V_0': 0.00036577 ,
            'Cp': [1862.2, -0.014018, -8983100.0, -14923.0] ,
            'a_0': 3.71e-05 ,
            'K_0': 51300000000.0 ,
            'Kprime_0': 7.33 ,
            'Kdprime_0': -1.43e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class atg (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg48.0Si34.0O147.0H62.0'
       self.params = {
            'name': 'atg',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -71404690.0 ,
            'S_0': 3620.0 ,
            'V_0': 0.0017548 ,
            'Cp': [9621.0, -0.091183, -35941600.0, -83034.2] ,
            'a_0': 2.8e-05 ,
            'K_0': 63100000000.0 ,
            'Kprime_0': 5.92 ,
            'Kdprime_0': -9.4e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class ab (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Na1.0Al1.0Si3.0O8.0'
       self.params = {
            'name': 'ab',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -3935480.0 ,
            'S_0': 207.4 ,
            'V_0': 0.00010067 ,
            'Cp': [452.0, -0.013364, -1275900.0, -3953.6] ,
            'a_0': 2.36e-05 ,
            'K_0': 54100000000.0 ,
            'Kprime_0': 5.91 ,
            'Kdprime_0': -1.09e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'BW_deltaH': 14000.0 ,
            'BW_deltaV': 4.2e-07 ,
            'BW_W': 13000.0 ,
            'BW_Wv': 4.2e-07 ,
            'BW_n': 3.0 ,
            'BW_factor': 0.9 }

class abh (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Na1.0Al1.0Si3.0O8.0'
       self.params = {
            'name': 'abh',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -3921480.0 ,
            'S_0': 224.3 ,
            'V_0': 0.00010105 ,
            'Cp': [452.0, -0.013364, -1275900.0, -3953.6] ,
            'a_0': 2.41e-05 ,
            'K_0': 54100000000.0 ,
            'Kprime_0': 5.91 ,
            'Kdprime_0': -1.09e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class mic (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='K1.0Al1.0Si3.0O8.0'
       self.params = {
            'name': 'mic',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -3975350.0 ,
            'S_0': 214.3 ,
            'V_0': 0.00010871 ,
            'Cp': [448.8, -0.010075, -1007300.0, -3973.1] ,
            'a_0': 1.66e-05 ,
            'K_0': 58300000000.0 ,
            'Kprime_0': 4.02 ,
            'Kdprime_0': -6.9e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class san (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='K1.0Al1.0Si3.0O8.0'
       self.params = {
            'name': 'san',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -3966700.0 ,
            'S_0': 214.3 ,
            'V_0': 0.00010871 ,
            'Cp': [448.8, -0.010075, -1007300.0, -3973.1] ,
            'a_0': 1.66e-05 ,
            'K_0': 58300000000.0 ,
            'Kprime_0': 4.02 ,
            'Kdprime_0': -6.9e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'BW_deltaH': 8650.0 ,
            'BW_deltaV': 2.4e-07 ,
            'BW_W': 8500.0 ,
            'BW_Wv': 2.4e-07 ,
            'BW_n': 3.0 ,
            'BW_factor': 0.8 }

class an (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca1.0Al2.0Si2.0O8.0'
       self.params = {
            'name': 'an',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -4232690.0 ,
            'S_0': 200.5 ,
            'V_0': 0.00010079 ,
            'Cp': [370.5, 0.01001, -4339100.0, -1960.6] ,
            'a_0': 1.41e-05 ,
            'K_0': 86000000000.0 ,
            'Kprime_0': 4.09 ,
            'Kdprime_0': -4.8e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'BW_deltaH': 42010.0 ,
            'BW_deltaV': 1e-06 ,
            'BW_W': 42000.0 ,
            'BW_Wv': 1e-06 ,
            'BW_n': 1.0 ,
            'BW_factor': 2.0 }

class kcm (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='K1.0Al1.0Si3.0O9.0H2.0'
       self.params = {
            'name': 'kcm',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -4232640.0 ,
            'S_0': 281.5 ,
            'V_0': 0.00011438 ,
            'Cp': [536.5, -0.01009, -980400.0, -4735.0] ,
            'a_0': 3.21e-05 ,
            'K_0': 42500000000.0 ,
            'Kprime_0': 2.0 ,
            'Kdprime_0': -4.7e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class wa (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='K2.0Si4.0O9.0'
       self.params = {
            'name': 'wa',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -4271890.0 ,
            'S_0': 254.0 ,
            'V_0': 0.00010844 ,
            'Cp': [499.1, 0.0, 0.0, -4350.1] ,
            'a_0': 2.66e-05 ,
            'K_0': 90000000000.0 ,
            'Kprime_0': 4.0 ,
            'Kdprime_0': -4.4e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class hol (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='K1.0Al1.0Si3.0O8.0'
       self.params = {
            'name': 'hol',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -3791960.0 ,
            'S_0': 166.2 ,
            'V_0': 7.128e-05 ,
            'Cp': [417.6, -0.003617, -4748100.0, -2819.9] ,
            'a_0': 2.8e-05 ,
            'K_0': 1.8e+11 ,
            'Kprime_0': 4.0 ,
            'Kdprime_0': -2.2e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class q (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Si1.0O2.0'
       self.params = {
            'name': 'q',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -910720.0 ,
            'S_0': 41.43 ,
            'V_0': 2.269e-05 ,
            'Cp': [92.9, -0.000642, -714900.0, -716.1] ,
            'a_0': 0.0 ,
            'K_0': 73000000000.0 ,
            'Kprime_0': 6.0 ,
            'Kdprime_0': -8.2e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'landau_Tc': 847.0 ,
            'landau_Smax': 4.95 ,
            'landau_Vmax': 1.188e-06 }

class trd (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Si1.0O2.0'
       self.params = {
            'name': 'trd',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -907110.0 ,
            'S_0': 44.1 ,
            'V_0': 2.8e-05 ,
            'Cp': [74.9, 0.0031, -1174000.0, -236.7] ,
            'a_0': 0.0 ,
            'K_0': 15000000000.0 ,
            'Kprime_0': 4.36 ,
            'Kdprime_0': -2.91e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class crst (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Si1.0O2.0'
       self.params = {
            'name': 'crst',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -904270.0 ,
            'S_0': 50.86 ,
            'V_0': 2.745e-05 ,
            'Cp': [72.7, 0.001304, -4129000.0, 0.0] ,
            'a_0': 0.0 ,
            'K_0': 16000000000.0 ,
            'Kprime_0': 4.35 ,
            'Kdprime_0': -2.72e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class coe (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Si1.0O2.0'
       self.params = {
            'name': 'coe',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -907000.0 ,
            'S_0': 39.6 ,
            'V_0': 2.064e-05 ,
            'Cp': [107.8, -0.003279, -190300.0, -1041.6] ,
            'a_0': 1.23e-05 ,
            'K_0': 97900000000.0 ,
            'Kprime_0': 4.19 ,
            'Kdprime_0': -4.3e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class stv (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Si1.0O2.0'
       self.params = {
            'name': 'stv',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -876390.0 ,
            'S_0': 24.0 ,
            'V_0': 1.401e-05 ,
            'Cp': [68.1, 0.00601, -1978200.0, -82.1] ,
            'a_0': 1.58e-05 ,
            'K_0': 3.09e+11 ,
            'Kprime_0': 4.6 ,
            'Kdprime_0': -1.5e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class ne (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Na1.0Al1.0Si1.0O4.0'
       self.params = {
            'name': 'ne',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2094560.0 ,
            'S_0': 124.4 ,
            'V_0': 5.419e-05 ,
            'Cp': [272.7, -0.012398, 0.0, -2763.1] ,
            'a_0': 4.63e-05 ,
            'K_0': 46500000000.0 ,
            'Kprime_0': 4.16 ,
            'Kdprime_0': -8.9e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'landau_Tc': 467.0 ,
            'landau_Smax': 10.0 ,
            'landau_Vmax': 8e-07 }

class cg (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Na1.0Al1.0Si1.0O4.0'
       self.params = {
            'name': 'cg',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2091720.0 ,
            'S_0': 118.7 ,
            'V_0': 5.603e-05 ,
            'Cp': [116.1, 0.086021, -1992700.0, 0.0] ,
            'a_0': 4.5e-05 ,
            'K_0': 46500000000.0 ,
            'Kprime_0': 4.16 ,
            'Kdprime_0': -8.9e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class cgh (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Na1.0Al1.0Si1.0O4.0'
       self.params = {
            'name': 'cgh',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2078010.0 ,
            'S_0': 135.0 ,
            'V_0': 5.67e-05 ,
            'Cp': [229.2, 0.011876, 0.0, -1970.7] ,
            'a_0': 4.67e-05 ,
            'K_0': 46500000000.0 ,
            'Kprime_0': 4.16 ,
            'Kdprime_0': -8.9e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class sdl (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Na8.0Al6.0Si6.0Cl2.0O24.0'
       self.params = {
            'name': 'sdl',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -13405530.0 ,
            'S_0': 910.0 ,
            'V_0': 0.0004213 ,
            'Cp': [1532.7, 0.047747, -2972800.0, -12427.0] ,
            'a_0': 4.63e-05 ,
            'K_0': 46500000000.0 ,
            'Kprime_0': 4.16 ,
            'Kdprime_0': -8.9e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class kls (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='K1.0Al1.0Si1.0O4.0'
       self.params = {
            'name': 'kls',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2122960.0 ,
            'S_0': 136.0 ,
            'V_0': 6.052e-05 ,
            'Cp': [242.0, -0.004482, -895800.0, -1935.8] ,
            'a_0': 3.16e-05 ,
            'K_0': 51400000000.0 ,
            'Kprime_0': 2.0 ,
            'Kdprime_0': -3.9e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class lc (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='K1.0Al1.0Si2.0O6.0'
       self.params = {
            'name': 'lc',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -3029270.0 ,
            'S_0': 198.5 ,
            'V_0': 8.826e-05 ,
            'Cp': [369.8, -0.016332, 684700.0, -3683.1] ,
            'a_0': 1.85e-05 ,
            'K_0': 45000000000.0 ,
            'Kprime_0': 5.7 ,
            'Kdprime_0': -1.27e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'BW_deltaH': 11610.0 ,
            'BW_deltaV': 4e-06 ,
            'BW_W': 11600.0 ,
            'BW_Wv': 4e-06 ,
            'BW_n': 2.0 ,
            'BW_factor': 0.7 }

class me (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca4.0Al6.0Si6.0O27.0C1.0'
       self.params = {
            'name': 'me',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -13841820.0 ,
            'S_0': 752.0 ,
            'V_0': 0.00033985 ,
            'Cp': [1359.0, 0.036442, -8594700.0, -9598.2] ,
            'a_0': 1.81e-05 ,
            'K_0': 87000000000.0 ,
            'Kprime_0': 4.09 ,
            'Kdprime_0': -4.7e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class wrk (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca1.0Al2.0Si4.0O14.0H4.0'
       self.params = {
            'name': 'wrk',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -6662450.0 ,
            'S_0': 380.0 ,
            'V_0': 0.0001904 ,
            'Cp': [838.3, -0.02146, -2272000.0, -7292.3] ,
            'a_0': 1.49e-05 ,
            'K_0': 86000000000.0 ,
            'Kprime_0': 4.09 ,
            'Kdprime_0': -4.8e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class lmt (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca1.0Al2.0Si4.0O16.0H8.0'
       self.params = {
            'name': 'lmt',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -7262700.0 ,
            'S_0': 465.0 ,
            'V_0': 0.0002037 ,
            'Cp': [1013.4, -0.021413, -2235800.0, -8806.7] ,
            'a_0': 1.37e-05 ,
            'K_0': 86000000000.0 ,
            'Kprime_0': 4.09 ,
            'Kdprime_0': -4.8e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class heu (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca1.0Al2.0Si7.0O24.0H12.0'
       self.params = {
            'name': 'heu',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -10545220.0 ,
            'S_0': 783.0 ,
            'V_0': 0.000317 ,
            'Cp': [1504.8, -0.033224, -2959300.0, -13297.2] ,
            'a_0': 1.57e-05 ,
            'K_0': 27400000000.0 ,
            'Kprime_0': 4.0 ,
            'Kdprime_0': -1.46e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class stlb (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca1.0Al2.0Si7.0O25.0H14.0'
       self.params = {
            'name': 'stlb',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -10896760.0 ,
            'S_0': 710.0 ,
            'V_0': 0.0003287 ,
            'Cp': [1588.4, -0.032043, -3071600.0, -13966.9] ,
            'a_0': 1.51e-05 ,
            'K_0': 86000000000.0 ,
            'Kprime_0': 4.09 ,
            'Kdprime_0': -4.8e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class anl (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Na1.0Al1.0Si2.0O7.0H2.0'
       self.params = {
            'name': 'anl',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -3307220.0 ,
            'S_0': 232.0 ,
            'V_0': 9.74e-05 ,
            'Cp': [643.5, -0.016067, 9302300.0, -9179.6] ,
            'a_0': 2.76e-05 ,
            'K_0': 40000000000.0 ,
            'Kprime_0': 4.18 ,
            'Kdprime_0': -1.04e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class lime (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca1.0O1.0'
       self.params = {
            'name': 'lime',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -634530.0 ,
            'S_0': 38.1 ,
            'V_0': 1.676e-05 ,
            'Cp': [52.4, 0.003673, -750700.0, -51.0] ,
            'a_0': 3.41e-05 ,
            'K_0': 1.13e+11 ,
            'Kprime_0': 3.87 ,
            'Kdprime_0': -3.4e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class ru (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ti1.0O2.0'
       self.params = {
            'name': 'ru',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -944360.0 ,
            'S_0': 50.5 ,
            'V_0': 1.882e-05 ,
            'Cp': [90.4, 0.0029, 0.0, -623.8] ,
            'a_0': 2.24e-05 ,
            'K_0': 2.22e+11 ,
            'Kprime_0': 4.24 ,
            'Kdprime_0': -1.9e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class per (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg1.0O1.0'
       self.params = {
            'name': 'per',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -601530.0 ,
            'S_0': 26.5 ,
            'V_0': 1.125e-05 ,
            'Cp': [60.5, 0.000362, -535800.0, -299.2] ,
            'a_0': 3.11e-05 ,
            'K_0': 1.616e+11 ,
            'Kprime_0': 3.95 ,
            'Kdprime_0': -2.4e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class fper (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe1.0O1.0'
       self.params = {
            'name': 'fper',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -259870.0 ,
            'S_0': 58.6 ,
            'V_0': 1.206e-05 ,
            'Cp': [44.4, 0.00828, -1214200.0, 185.2] ,
            'a_0': 3.22e-05 ,
            'K_0': 1.52e+11 ,
            'Kprime_0': 4.9 ,
            'Kdprime_0': -3.2e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class mang (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mn1.0O1.0'
       self.params = {
            'name': 'mang',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -385550.0 ,
            'S_0': 59.7 ,
            'V_0': 1.322e-05 ,
            'Cp': [59.8, 0.0036, -31400.0, -282.6] ,
            'a_0': 3.69e-05 ,
            'K_0': 1.645e+11 ,
            'Kprime_0': 4.46 ,
            'Kdprime_0': -2.7e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class cor (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Al2.0O3.0'
       self.params = {
            'name': 'cor',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1675270.0 ,
            'S_0': 50.9 ,
            'V_0': 2.558e-05 ,
            'Cp': [139.5, 0.00589, -2460600.0, -589.2] ,
            'a_0': 1.8e-05 ,
            'K_0': 2.54e+11 ,
            'Kprime_0': 4.34 ,
            'Kdprime_0': -1.7e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class mcor (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg1.0Si1.0O3.0'
       self.params = {
            'name': 'mcor',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1474440.0 ,
            'S_0': 59.3 ,
            'V_0': 2.635e-05 ,
            'Cp': [147.8, 0.002015, -2395000.0, -801.8] ,
            'a_0': 2.12e-05 ,
            'K_0': 2.11e+11 ,
            'Kprime_0': 4.55 ,
            'Kdprime_0': -2.2e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class hem (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe2.0O3.0'
       self.params = {
            'name': 'hem',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -825610.0 ,
            'S_0': 87.4 ,
            'V_0': 3.027e-05 ,
            'Cp': [163.9, 0.0, -2257200.0, -657.6] ,
            'a_0': 2.79e-05 ,
            'K_0': 2.23e+11 ,
            'Kprime_0': 4.04 ,
            'Kdprime_0': -1.8e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'landau_Tc': 955.0 ,
            'landau_Smax': 15.6 ,
            'landau_Vmax': 0.0 }

class esk (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Cr2.0O3.0'
       self.params = {
            'name': 'esk',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1137320.0 ,
            'S_0': 83.0 ,
            'V_0': 2.909e-05 ,
            'Cp': [119.0, 0.009496, -1442000.0, -3.4] ,
            'a_0': 1.59e-05 ,
            'K_0': 2.38e+11 ,
            'Kprime_0': 4.0 ,
            'Kdprime_0': -1.7e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class bix (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mn2.0O3.0'
       self.params = {
            'name': 'bix',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -959000.0 ,
            'S_0': 113.7 ,
            'V_0': 3.137e-05 ,
            'Cp': [145.1, 0.023534, 721600.0, -1008.4] ,
            'a_0': 2.91e-05 ,
            'K_0': 2.23e+11 ,
            'Kprime_0': 4.04 ,
            'Kdprime_0': -1.8e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class NiO (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ni1.0O1.0'
       self.params = {
            'name': 'NiO',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -239470.0 ,
            'S_0': 38.0 ,
            'V_0': 1.097e-05 ,
            'Cp': [47.7, 0.007824, -392500.0, 0.0] ,
            'a_0': 3.3e-05 ,
            'K_0': 2e+11 ,
            'Kprime_0': 3.94 ,
            'Kdprime_0': -2e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'landau_Tc': 520.0 ,
            'landau_Smax': 5.7 ,
            'landau_Vmax': 0.0 }

class pnt (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mn1.0Ti1.0O3.0'
       self.params = {
            'name': 'pnt',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1361950.0 ,
            'S_0': 105.5 ,
            'V_0': 3.288e-05 ,
            'Cp': [143.5, 0.003373, -1940700.0, -407.6] ,
            'a_0': 2.4e-05 ,
            'K_0': 1.7e+11 ,
            'Kprime_0': 8.3 ,
            'Kdprime_0': -4.9e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class geik (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg1.0Ti1.0O3.0'
       self.params = {
            'name': 'geik',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1568960.0 ,
            'S_0': 73.6 ,
            'V_0': 3.086e-05 ,
            'Cp': [151.0, 0.0, -1890400.0, -652.2] ,
            'a_0': 2.15e-05 ,
            'K_0': 1.7e+11 ,
            'Kprime_0': 8.3 ,
            'Kdprime_0': -4.9e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class ilm (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe1.0Ti1.0O3.0'
       self.params = {
            'name': 'ilm',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1230450.0 ,
            'S_0': 109.5 ,
            'V_0': 3.169e-05 ,
            'Cp': [138.9, 0.005081, -1288800.0, -463.7] ,
            'a_0': 2.4e-05 ,
            'K_0': 1.7e+11 ,
            'Kprime_0': 8.3 ,
            'Kdprime_0': -4.9e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'landau_Tc': 1900.0 ,
            'landau_Smax': 12.0 ,
            'landau_Vmax': 2e-07 }

class bdy (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Zr1.0O2.0'
       self.params = {
            'name': 'bdy',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1100340.0 ,
            'S_0': 50.4 ,
            'V_0': 2.115e-05 ,
            'Cp': [103.5, -0.004547, -416200.0, -713.6] ,
            'a_0': 2e-05 ,
            'K_0': 95300000000.0 ,
            'Kprime_0': 3.88 ,
            'Kdprime_0': -4.1e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class ten (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Cu1.0O1.0'
       self.params = {
            'name': 'ten',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -156100.0 ,
            'S_0': 42.6 ,
            'V_0': 1.222e-05 ,
            'Cp': [31.0, 0.01374, -1258000.0, 369.3] ,
            'a_0': 3.57e-05 ,
            'K_0': 2e+11 ,
            'Kprime_0': 3.94 ,
            'Kdprime_0': -2e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class cup (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Cu2.0O1.0'
       self.params = {
            'name': 'cup',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -170600.0 ,
            'S_0': 92.4 ,
            'V_0': 2.344e-05 ,
            'Cp': [110.3, 0.0, 0.0, -674.8] ,
            'a_0': 3.33e-05 ,
            'K_0': 1.31e+11 ,
            'Kprime_0': 5.7 ,
            'Kdprime_0': -4.3e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class sp (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg1.0Al2.0O4.0'
       self.params = {
            'name': 'sp',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2301190.0 ,
            'S_0': 82.0 ,
            'V_0': 3.978e-05 ,
            'Cp': [222.9, 0.006127, -1686000.0, -1551.0] ,
            'a_0': 1.93e-05 ,
            'K_0': 1.922e+11 ,
            'Kprime_0': 4.04 ,
            'Kdprime_0': -2.1e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'BW_deltaH': 8000.0 ,
            'BW_deltaV': 0.0 ,
            'BW_W': 1200.0 ,
            'BW_Wv': 0.0 ,
            'BW_n': 2.0 ,
            'BW_factor': 0.5 }

class herc (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe1.0Al2.0O4.0'
       self.params = {
            'name': 'herc',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1953030.0 ,
            'S_0': 113.9 ,
            'V_0': 4.075e-05 ,
            'Cp': [216.7, 0.005868, -2430200.0, -1178.3] ,
            'a_0': 2.06e-05 ,
            'K_0': 1.922e+11 ,
            'Kprime_0': 4.04 ,
            'Kdprime_0': -2.1e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'BW_deltaH': 18300.0 ,
            'BW_deltaV': 0.0 ,
            'BW_W': 13600.0 ,
            'BW_Wv': 0.0 ,
            'BW_n': 2.0 ,
            'BW_factor': 1.0 }

class mt (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe3.0O4.0'
       self.params = {
            'name': 'mt',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1114500.0 ,
            'S_0': 146.9 ,
            'V_0': 4.452e-05 ,
            'Cp': [262.5, -0.007205, -1926200.0, -1655.7] ,
            'a_0': 3.71e-05 ,
            'K_0': 1.857e+11 ,
            'Kprime_0': 4.05 ,
            'Kdprime_0': -2.2e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'landau_Tc': 848.0 ,
            'landau_Smax': 35.0 ,
            'landau_Vmax': 0.0 }

class mft (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg1.0Fe2.0O4.0'
       self.params = {
            'name': 'mft',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1442290.0 ,
            'S_0': 121.0 ,
            'V_0': 4.457e-05 ,
            'Cp': [270.5, -0.007505, -999200.0, -2022.4] ,
            'a_0': 3.63e-05 ,
            'K_0': 1.857e+11 ,
            'Kprime_0': 4.05 ,
            'Kdprime_0': -2.2e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'landau_Tc': 665.0 ,
            'landau_Smax': 17.0 ,
            'landau_Vmax': 0.0 }

class usp (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe2.0Ti1.0O4.0'
       self.params = {
            'name': 'usp',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1491120.0 ,
            'S_0': 180.0 ,
            'V_0': 4.682e-05 ,
            'Cp': [-102.6, 0.14252, -9144500.0, 5270.7] ,
            'a_0': 3.86e-05 ,
            'K_0': 1.857e+11 ,
            'Kprime_0': 4.05 ,
            'Kdprime_0': -2.2e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class picr (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg1.0Cr2.0O4.0'
       self.params = {
            'name': 'picr',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1762600.0 ,
            'S_0': 118.3 ,
            'V_0': 4.356e-05 ,
            'Cp': [196.1, 0.005398, -3126000.0, -616.9] ,
            'a_0': 1.8e-05 ,
            'K_0': 1.922e+11 ,
            'Kprime_0': 4.04 ,
            'Kdprime_0': -2.1e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'BW_deltaH': 8000.0 ,
            'BW_deltaV': 0.0 ,
            'BW_W': 1200.0 ,
            'BW_Wv': 0.0 ,
            'BW_n': 2.0 ,
            'BW_factor': 0.5 }

class br (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg1.0O2.0H2.0'
       self.params = {
            'name': 'br',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -925560.0 ,
            'S_0': 63.2 ,
            'V_0': 2.463e-05 ,
            'Cp': [158.4, -0.004076, -1052300.0, -1171.3] ,
            'a_0': 6.2e-05 ,
            'K_0': 41500000000.0 ,
            'Kprime_0': 6.45 ,
            'Kdprime_0': -1.55e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class dsp (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Al1.0O2.0H1.0'
       self.params = {
            'name': 'dsp',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -999840.0 ,
            'S_0': 34.5 ,
            'V_0': 1.786e-05 ,
            'Cp': [145.1, 0.008709, 584400.0, -1741.1] ,
            'a_0': 3.57e-05 ,
            'K_0': 2.28e+11 ,
            'Kprime_0': 4.04 ,
            'Kdprime_0': -1.8e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class gth (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe1.0O2.0H1.0'
       self.params = {
            'name': 'gth',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -561770.0 ,
            'S_0': 60.3 ,
            'V_0': 2.082e-05 ,
            'Cp': [139.3, 0.000147, -212700.0, -1077.8] ,
            'a_0': 4.35e-05 ,
            'K_0': 2.5e+11 ,
            'Kprime_0': 4.03 ,
            'Kdprime_0': -1.6e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class cc (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca1.0C1.0O3.0'
       self.params = {
            'name': 'cc',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1207760.0 ,
            'S_0': 92.5 ,
            'V_0': 3.689e-05 ,
            'Cp': [140.9, 0.005029, -950700.0, -858.4] ,
            'a_0': 2.52e-05 ,
            'K_0': 73300000000.0 ,
            'Kprime_0': 4.06 ,
            'Kdprime_0': -5.5e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'landau_Tc': 1240.0 ,
            'landau_Smax': 10.0 ,
            'landau_Vmax': 4e-07 }

class arag (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca1.0O3.0C1.0'
       self.params = {
            'name': 'arag',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1207650.0 ,
            'S_0': 89.8 ,
            'V_0': 3.415e-05 ,
            'Cp': [167.1, 0.010695, 162000.0, -1564.9] ,
            'a_0': 6.14e-05 ,
            'K_0': 61400000000.0 ,
            'Kprime_0': 5.87 ,
            'Kdprime_0': -9.6e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class mag (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mg1.0O3.0C1.0'
       self.params = {
            'name': 'mag',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1110920.0 ,
            'S_0': 65.5 ,
            'V_0': 2.803e-05 ,
            'Cp': [186.4, -0.003772, 0.0, -1886.2] ,
            'a_0': 3.38e-05 ,
            'K_0': 1.028e+11 ,
            'Kprime_0': 5.41 ,
            'Kdprime_0': -5.3e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class sid (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe1.0C1.0O3.0'
       self.params = {
            'name': 'sid',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -762220.0 ,
            'S_0': 93.3 ,
            'V_0': 2.943e-05 ,
            'Cp': [168.4, 0.0, 0.0, -1483.6] ,
            'a_0': 4.39e-05 ,
            'K_0': 1.2e+11 ,
            'Kprime_0': 4.07 ,
            'Kdprime_0': -3.4e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class rhc (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Mn1.0C1.0O3.0'
       self.params = {
            'name': 'rhc',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -892280.0 ,
            'S_0': 98.0 ,
            'V_0': 3.107e-05 ,
            'Cp': [169.5, 0.0, 0.0, -1534.3] ,
            'a_0': 2.44e-05 ,
            'K_0': 95300000000.0 ,
            'Kprime_0': 3.88 ,
            'Kdprime_0': -4.1e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class dol (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca1.0Mg1.0O6.0C2.0'
       self.params = {
            'name': 'dol',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2326220.0 ,
            'S_0': 156.1 ,
            'V_0': 6.429e-05 ,
            'Cp': [358.9, -0.004905, 0.0, -3456.2] ,
            'a_0': 3.28e-05 ,
            'K_0': 94300000000.0 ,
            'Kprime_0': 3.74 ,
            'Kdprime_0': -4e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'BW_deltaH': 11910.0 ,
            'BW_deltaV': 1.6e-07 ,
            'BW_W': 11900.0 ,
            'BW_Wv': 1.6e-07 ,
            'BW_n': 1.0 ,
            'BW_factor': 1.0 }

class ank (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca1.0Fe1.0O6.0C2.0'
       self.params = {
            'name': 'ank',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1971410.0 ,
            'S_0': 188.46 ,
            'V_0': 6.606e-05 ,
            'Cp': [341.0, -0.001161, 0.0, -3054.8] ,
            'a_0': 3.46e-05 ,
            'K_0': 91400000000.0 ,
            'Kprime_0': 3.88 ,
            'Kdprime_0': -4.3e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'BW_deltaH': 11910.0 ,
            'BW_deltaV': 1.6e-07 ,
            'BW_W': 11900.0 ,
            'BW_Wv': 1.6e-07 ,
            'BW_n': 1.0 ,
            'BW_factor': 1.0 }

class syv (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='K1.0Cl1.0'
       self.params = {
            'name': 'syv',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -436500.0 ,
            'S_0': 82.6 ,
            'V_0': 3.752e-05 ,
            'Cp': [46.2, 0.01797, 0.0, 0.0] ,
            'a_0': 0.0001109 ,
            'K_0': 17000000000.0 ,
            'Kprime_0': 5.0 ,
            'Kdprime_0': -2.94e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class hlt (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Na1.0Cl1.0'
       self.params = {
            'name': 'hlt',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -411300.0 ,
            'S_0': 72.1 ,
            'V_0': 2.702e-05 ,
            'Cp': [45.2, 0.01797, 0.0, 0.0] ,
            'a_0': 0.0001147 ,
            'K_0': 23800000000.0 ,
            'Kprime_0': 5.0 ,
            'Kdprime_0': -2.1e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class pyr (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe1.0S2.0'
       self.params = {
            'name': 'pyr',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -171640.0 ,
            'S_0': 52.9 ,
            'V_0': 2.394e-05 ,
            'Cp': [37.3, 0.026715, -1817000.0, 649.3] ,
            'a_0': 3.1e-05 ,
            'K_0': 1.395e+11 ,
            'Kprime_0': 4.09 ,
            'Kdprime_0': -2.9e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class trot (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe1.0S1.0'
       self.params = {
            'name': 'trot',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -99030.0 ,
            'S_0': 65.5 ,
            'V_0': 1.819e-05 ,
            'Cp': [50.2, 0.011052, -940000.0, 0.0] ,
            'a_0': 5.68e-05 ,
            'K_0': 65800000000.0 ,
            'Kprime_0': 4.17 ,
            'Kdprime_0': -6.3e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'landau_Tc': 598.0 ,
            'landau_Smax': 12.0 ,
            'landau_Vmax': 4.1e-07 }

class tro (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe1.0S1.0'
       self.params = {
            'name': 'tro',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -97760.0 ,
            'S_0': 70.8 ,
            'V_0': 1.819e-05 ,
            'Cp': [50.2, 0.011052, -940000.0, 0.0] ,
            'a_0': 5.73e-05 ,
            'K_0': 65800000000.0 ,
            'Kprime_0': 4.17 ,
            'Kdprime_0': -6.3e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'landau_Tc': 598.0 ,
            'landau_Smax': 12.0 ,
            'landau_Vmax': 4.1e-07 }

class lot (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe1.0S1.0'
       self.params = {
            'name': 'lot',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -102160.0 ,
            'S_0': 60.0 ,
            'V_0': 1.818e-05 ,
            'Cp': [50.2, 0.011052, -940000.0, 0.0] ,
            'a_0': 4.93e-05 ,
            'K_0': 65800000000.0 ,
            'Kprime_0': 4.17 ,
            'Kdprime_0': -6.3e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'landau_Tc': 420.0 ,
            'landau_Smax': 10.0 ,
            'landau_Vmax': 0.0 }

class trov (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe0.875S1.0'
       self.params = {
            'name': 'trov',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -96020.0 ,
            'S_0': 57.5 ,
            'V_0': 1.738e-05 ,
            'Cp': [51.1, 0.008307, -669700.0, 0.0] ,
            'a_0': 5.94e-05 ,
            'K_0': 65800000000.0 ,
            'Kprime_0': 4.17 ,
            'Kdprime_0': -6.3e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'landau_Tc': 595.0 ,
            'landau_Smax': 10.0 ,
            'landau_Vmax': 1.6e-07 }

class any (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ca1.0S1.0O4.0'
       self.params = {
            'name': 'any',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1434400.0 ,
            'S_0': 106.9 ,
            'V_0': 4.594e-05 ,
            'Cp': [128.7, 0.048545, -1223000.0, -560.5] ,
            'a_0': 4.18e-05 ,
            'K_0': 54380000000.0 ,
            'Kprime_0': 4.19 ,
            'Kdprime_0': -7.7e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class iron (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Fe1.0'
       self.params = {
            'name': 'iron',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -0.0 ,
            'S_0': 27.09 ,
            'V_0': 7.09e-06 ,
            'Cp': [46.2, 0.005159, 723100.0, -556.2] ,
            'a_0': 3.56e-05 ,
            'K_0': 1.64e+11 ,
            'Kprime_0': 5.16 ,
            'Kdprime_0': -3.1e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'landau_Tc': 1042.0 ,
            'landau_Smax': 8.3 ,
            'landau_Vmax': 0.0 }

class Ni (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Ni1.0'
       self.params = {
            'name': 'Ni',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': 0.0 ,
            'S_0': 29.87 ,
            'V_0': 6.59e-06 ,
            'Cp': [49.8, 0.0, 585900.0, -533.9] ,
            'a_0': 4.28e-05 ,
            'K_0': 1.905e+11 ,
            'Kprime_0': 4.25 ,
            'Kdprime_0': -2.2e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'landau_Tc': 631.0 ,
            'landau_Smax': 3.0 ,
            'landau_Vmax': 0.0 }

class Cu (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='Cu1.0'
       self.params = {
            'name': 'Cu',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -0.0 ,
            'S_0': 33.14 ,
            'V_0': 7.11e-06 ,
            'Cp': [12.4, 0.00922, -379900.0, 233.5] ,
            'a_0': 3.58e-05 ,
            'K_0': 1.625e+11 ,
            'Kprime_0': 4.24 ,
            'Kdprime_0': -2.6e-11 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class gph (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='C1.0'
       self.params = {
            'name': 'gph',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': 0.0 ,
            'S_0': 5.76 ,
            'V_0': 5.3e-06 ,
            'Cp': [34.3, 0.0, -240700.0, -403.8] ,
            'a_0': 1.65e-05 ,
            'K_0': 31200000000.0 ,
            'Kprime_0': 3.9 ,
            'Kdprime_0': -1.25e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class diam (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='C1.0'
       self.params = {
            'name': 'diam',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': 1890.0 ,
            'S_0': 2.36 ,
            'V_0': 3.42e-06 ,
            'Cp': [40.0, 0.0, -28500.0, -580.5] ,
            'a_0': 4e-06 ,
            'K_0': 4.465e+11 ,
            'Kprime_0': 1.61 ,
            'Kdprime_0': -3.6e-12 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

class S (Mineral):
    """
    Holland and Powell, 2011, and references therein
    """
    def __init__(self):
       formula='S1.0'
       self.params = {
            'name': 'S',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': 0.0 ,
            'S_0': 32.05 ,
            'V_0': 1.551e-05 ,
            'Cp': [56.6, -0.004557, 638000.0, -681.8] ,
            'a_0': 6.4e-05 ,
            'K_0': 14500000000.0 ,
            'Kprime_0': 7.0 ,
            'Kdprime_0': -4.8e-10 ,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}


