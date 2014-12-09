# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import burnman.mineral_helpers as bmb
from burnman.mineral import Mineral
import numpy as np
                
class forsterite (Mineral):
    """
        mg olivine
        Stixrude & Lithgow-Bertelloni 2011 and references therein
        """
    def __init__(self,perturb=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state': 'slb3',
            'V_0': 43.60e-6,
            'K_0': 128.0e9,
            'err_K_0':2.e9,
            'Kprime_0': 4.2,
            'err_Kprime_0':0.2,
            'G_0': 82.0e9,
            'err_G_0':2.e9,
            'Gprime_0': 1.5,
            'err_Gprime_0':0.1,
            'molar_mass': .14069,
            'n': 7,
            'Debye_0': 809.,
            'err_Debye_0' : 1.,
            'grueneisen_0': .99,
            'err_grueneisen_0': .03,
            'q_0': 2.1,
            'err_q_0': 0.2,
            'eta_s_0': 2.3,
            'err_eta_s_0' : 0.1 }


class fayalite (Mineral):
    """
        fe olivine
        Stixrude & Lithgow-Bertelloni 2011 and references therein
        """
    def __init__(self,perturb=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state': 'slb3',
            'V_0': 46.29e-6,
            'K_0': 135.0e9,
            'err_K_0':2e9,
            'Kprime_0': 4.2,
            'err_Kprime_0':1.0,
            'G_0': 51.0e9,
            'err_G_0':2.e9,
            'Gprime_0': 1.5,
            'err_Gprime_0':0.5,
            'molar_mass': .20377,
            'n': 7,
            'Debye_0': 619.,
            'err_Debye_0' : 2.,
            'grueneisen_0': 1.06,
            'err_grueneisen_0': .07,
            'q_0': 3.6,
            'err_q_0': 1.0,
            'eta_s_0': 1.0,
            'err_eta_s_0' : 0.6 }

class mg_wadsleyite (Mineral):
    """
        Stixrude & Lithgow-Bertelloni 2011 and references therein
        """
    def __init__(self,perturb=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state': 'slb3',
            'V_0': 40.52e-6,
            'K_0': 169.0e9,
            'err_K_0':3.e9,
            'Kprime_0': 4.3,
            'err_Kprime_0':0.2,
            'G_0': 112.0e9,
            'err_G_0':2.e9,
            'Gprime_0': 1.4, #1.4
            'err_Gprime_0':0.2,
            'molar_mass': .14069,
            'n': 7,
            'Debye_0': 844.,
            'err_Debye_0' : 7.,
            'grueneisen_0': 1.21,
            'err_grueneisen_0': .09,
            'q_0': 2.0,
            'err_q_0': 1.0,
            'eta_s_0': 2.6,
            'err_eta_s_0' : 0.4 }


class fe_wadsleyite (Mineral):
    """
        Stixrude & Lithgow-Bertelloni 2011 and references therein
        """
    def __init__(self,perturb=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state': 'slb3',
            'V_0': 42.8e-6,
            'K_0': 169.0e9,
            'err_K_0':13.e9,
            'Kprime_0': 4.3,
            'err_Kprime_0':1.0,
            'G_0': 72.0e9, #72.e9
            'err_G_0':12.e9,
            'Gprime_0': 1.4, #1.4
            'err_Gprime_0':0.5,
            'molar_mass': .20377,
            'n': 7,
            'Debye_0': 665.,
            'err_Debye_0' : 21.,
            'grueneisen_0': 1.21,
            'err_grueneisen_0': .3,
            'q_0': 2.0,
            'err_q_0': 1.0,
            'eta_s_0': 1.0,
            'err_eta_s_0' : 1.0 }


class mg_ringwoodite (Mineral):
    """
        Stixrude & Lithgow-Bertelloni 2011 and references therein
        """
    def __init__(self,perturb=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state': 'slb3',
            'V_0': 39.49e-6,
            'K_0': 185.0e9,
            'err_K_0':2.e9,
            'Kprime_0': 4.2,
            'err_Kprime_0':0.2,
            'G_0': 123.0e9,
            'err_G_0':2.e9,
            'Gprime_0': 1.4, #1.4
            'err_Gprime_0':0.1,
            'molar_mass': .14069,
            'n': 7,
            'Debye_0': 878.,
            'err_Debye_0' : 8.,
            'grueneisen_0': 1.11,
            'err_grueneisen_0': .1,
            'q_0': 2.4,
            'err_q_0': .4,
            'eta_s_0': 2.3,
            'err_eta_s_0' : 0.5 }


class fe_ringwoodite (Mineral):
    """
        Stixrude & Lithgow-Bertelloni 2011 and references therein
        """
    def __init__(self,perturb=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state': 'slb3',
            'V_0': 41.86e-6,
            'K_0': 213.0e9,
            'err_K_0':7.e9,
            'Kprime_0': 4.2,
            'err_Kprime_0':1.0,
            'G_0': 92.0e9, #92.e9
            'err_G_0':10.e9,
            'Gprime_0': 1.4,#1.4
            'err_Gprime_0':0.5,
            'molar_mass': .20377,
            'n': 7,
            'Debye_0': 679.,
            'err_Debye_0' : 8.,
            'grueneisen_0': 1.27,
            'err_grueneisen_0': .23,
            'q_0': 2.4,
            'err_q_0': 1.0,
            'eta_s_0': 1.8,
            'err_eta_s_0' : 1.0 }

'''
class mg_olivine_MTZ_ugly1(helper_phase_transitions):
    """
        stitches together differen MG_endmember components through the mantle transition zone. This is extremely simplified. Clapeyron slopes for the 410, 520, and 660, respectively, are taken from Katsuro et al. 2004, Inoue et al. 2006, and Fei et al. 2004
    """
    def __init__(self,fe_num,perturb=np.zeros([4,9])):
        self.perturb=perturb
        self.fe_num=fe_num

    def set_phase(self,pressure, temperature):

        #if ( pressure >= 25.47e9-1.3e6*temperature): #Fei et al. 2004
        if ( pressure >= 28.71e9-2.9e6*temperature): #Yu et al. 2007
            mat = mg_fe_perovskite(self.fe_num)
        else:
            if ( pressure >= 13.1e9+4.11e6*temperature ): #Inoue et al. 2006
                mat = mg_fe_ringwoodite(self.fe_num)
            else:
                if ( pressure >= 7.8e9+4.e6*temperature) :   #Katsuro et al. 2004
                    mat = mg_fe_wadsleyite(self.fe_num)
                else:
                    mat = mg_fe_olivine(self.fe_num)
        self.params=mat.params


class mg_olivine_MTZ_ugly2(helper_phase_transitions):
    """
        stitches together differen MG_endmember components through the mantle transition zone. This is extremely simplified. Clapeyron slopes for the 410, 520, and 660, respectively, are taken from Katsuro et al. 2004, Inoue et al. 2006, and Fei et al. 2004
        """
    def __init__(self,fe_num,perturb=np.zeros([4,9])):
        self.perturb=perturb
        self.fe_num=fe_num
    
    def set_phase(self,pressure, temperature,):

        #if ( pressure >= 25.47e9-1.3e6*temperature): #Fei et al. 2004
        if ( pressure >= 28.71e9-2.9e6*temperature): #Yu et al. 2007
            mat = ferropericlase(self.fe_num)
        else:
            if ( pressure >= 13.1e9+4.11e6*temperature ):
                mat = mg_fe_ringwoodite(self.fe_num)
            else:
                if ( pressure >= 7.8e9+4.e6*temperature) :
                    mat = mg_fe_wadsleyite(self.fe_num)
                else:
                    mat = mg_fe_olivine(self.fe_num)
        self.params=mat.params
'''




class stishovite (Mineral):
    """
    Stixrude & Lithgow-Bertelloni 2011 and references therein 
    """
    def __init__(self,perturb=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state': 'slb3',
            'V_0': 14.02e-6,
            'K_0': 314.0e9,
            'err_K_0':8.e9,
            'Kprime_0': 3.8,
            'err_Kprime_0':0.1,
            'G_0': 220.0e9,
            'err_G_0':12.e9,
            'Gprime_0': 1.9,
            'err_Gprime_0':0.1,
            'molar_mass': .0601,
            'n': 3,
            'Debye_0': 1108.,
            'err_Debye_0' : 13.,
            'grueneisen_0': 1.37,
            'err_grueneisen_0': .17,
            'q_0': 2.8,
            'err_q_0': 2.2,    # decreased so things don't crash... not the published value (which is 2.2)
            'eta_s_0': 4.6,
            'err_eta_s_0' : 1.0 }



class periclase (Mineral):
    """
    Stixrude & Lithgow-Bertelloni 2011 and references therein 
    """
    def __init__(self,perturb=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 11.24e-6,
            'K_0': 161.0e9,
            'err_K_0': 3.e9,
            'Kprime_0': 3.8,
            'err_Kprime_0':.2,
            'G_0': 131.0e9,
            'err_G_0':1.0e9,
            'Gprime_0': 2.1,
            'err_Gprime_0':.1,
            'molar_mass': .0403,
            'n': 2,
            'Debye_0': 767.,
            'err_Debye_0':9.,
            'grueneisen_0': 1.36,
            'err_grueneisen_0':.05,
            'q_0': 1.7,
            'err_q_0':.2,
            'eta_s_0': 2.8,
            'err_eta_s_0':.2 }


class wuestite (Mineral):
    """
    Stixrude & Lithgow-Bertelloni 2011 and references therein 
    """
    def __init__(self,perturb=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 12.26e-6,
            'K_0': 179.0e9,
            'err_K_0':1.e9,
            'Kprime_0': 4.9,
            'err_Kprime_0':.2,
            'G_0': 59.0e9,
            'err_G_0':1.e9,
            'Gprime_0': 1.4,
            'err_Gprime_0':.1,
            'molar_mass': .0718,
            'n': 2,
            'Debye_0': 454.,
            'err_Debye_0':21.,
            'grueneisen_0': 1.53,
            'err_grueneisen_0':.13,
            'q_0': 1.7,
            'err_q_0':1.0,
            'eta_s_0': -0.1,
            'err_eta_s_0':1.0}



class ferropericlase(bmb.HelperSolidSolution):
    def __init__(self, fe_num):
        base_materials = [periclase(), wuestite()]
        molar_fraction = [1. - fe_num, 0.0 + fe_num] # keep the 0.0 +, otherwise it is an array sometimes
        bmb.HelperSolidSolution.__init__(self, base_materials, molar_fraction)


class mg_fe_perovskite(bmb.HelperSolidSolution):
    def __init__(self, fe_num):
        base_materials = [mg_perovskite(), fe_perovskite()]
        molar_fraction = [1. - fe_num, 0.0 + fe_num] # keep the 0.0 +, otherwise it is an array sometimes
        bmb.HelperSolidSolution.__init__(self, base_materials, molar_fraction)

class mg_fe_olivine(bmb.HelperSolidSolution):
    def __init__(self, fe_num):
        base_materials = [forsterite(), fayalite()]
        molar_fraction = [1. - fe_num, 0.0 + fe_num] # keep the 0.0 +, otherwise it is an array sometimes
        bmb.HelperSolidSolution.__init__(self, base_materials, molar_fraction)

class mg_fe_wadsleyite(bmb.HelperSolidSolution):
    def __init__(self, fe_num):
        base_materials = [mg_wadsleyite(), fe_wadsleyite()]
        molar_fraction = [1. - fe_num, 0.0 + fe_num] # keep the 0.0 +, otherwise it is an array sometimes
        bmb.HelperSolidSolution.__init__(self, base_materials, molar_fraction)

class mg_fe_ringwoodite(bmb.HelperSolidSolution):
    def __init__(self, fe_num):
        base_materials = [mg_ringwoodite(), fe_ringwoodite()]
        molar_fraction = [1. - fe_num, 0.0 + fe_num] # keep the 0.0 +, otherwise it is an array sometimes
        bmb.HelperSolidSolution.__init__(self, base_materials, molar_fraction)

class mg_perovskite(Mineral):
    """
    Stixrude & Lithgow-Bertelloni 2011 and references therein  
    """
    def __init__(self,perturb=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 24.45e-6,
            'K_0': 251.0e9,
            'err_K_0': 3.e9,
            'Kprime_0': 4.1,
            'err_Kprime_0': 0.1,
            'G_0': 173.0e9,
            'err_G_0': 2.e9,
            'Gprime_0': 1.7,
            'err_Gprime_0' : 0.0,
            'molar_mass': .1000,
            'n': 5,
            'Debye_0': 905.,
            'err_Debye_0': 5.,
            'grueneisen_0': 1.57,
            'err_grueneisen_0':.05,
            'q_0': 1.1,
            'err_q_0': .3,
            'eta_s_0': 2.6,
            'err_eta_s_0':.3}


class fe_perovskite(Mineral):
    """
    Stixrude & Lithgow-Bertelloni 2011 and references therein 
    """
    def __init__(self,perturb=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 25.49e-6,
            'K_0': 272.0e9,
            'err_K_0':40e9,
            'Kprime_0': 4.1,
            'err_Kprime_0':1.,
            'G_0': 133.0e9,
            'err_G_0':40e9,
            'Gprime_0': 1.4,
            'err_Gprime_0':0.0,
            'molar_mass': .1319, 
            'n': 5,
            'Debye_0': 871.,
            'err_Debye_0':26.,
            'grueneisen_0': 1.57,
            'err_grueneisen_0':.3,
            'q_0': 1.1,
            'err_q_0':1.0,
            'eta_s_0': 2.3,
            'err_eta_s_0':1.0}




class mg_fe_perovskite_pt_dependent(bmb.HelperFeDependent):
    def __init__(self, iron_number_with_pt, idx):
        bmb.HelperFeDependent.__init__(self, iron_number_with_pt, idx)

    def create_inner_material(self, iron_number):
        return mg_fe_perovskite(iron_number)

class ferropericlase_pt_dependent(bmb.HelperFeDependent):
    def __init__(self, iron_number_with_pt, idx):
        bmb.HelperFeDependent.__init__(self, iron_number_with_pt, idx)

    def create_inner_material(self, iron_number):
        return ferropericlase(iron_number)


class pyrope(Mineral):
    """
        garnet phase
        Stixrude & Lithgow-Bertelloni 2011 and references therein
        """
    def __init__(self,perturb=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 113.08e-6,
            'K_0': 170.0e9,
            'err_K_0':2e9,
            'Kprime_0': 4.1,
            'err_Kprime_0':.3,
            'G_0': 94.0e9,
            'err_G_0':2e9,
            'Gprime_0': 1.4,
            'err_Gprime_0':0.2,
            'molar_mass': .4031,
            'n': 20,
            'Debye_0': 823.,
            'err_Debye_0':4.,
            'grueneisen_0': 1.01,
            'err_grueneisen_0':.06,
            'q_0': 1.4,
            'err_q_0':0.5,
            'eta_s_0': 1.0,
            'err_eta_s_0':0.3}

class grossular(Mineral):
    """
        garnet phase
        Stixrude & Lithgow-Bertelloni 2011 and references therein
        """
    def __init__(self,perturb=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 125.12e-6,
            'K_0': 167.0e9,
            'err_K_0':1e9,
            'Kprime_0': 3.9,
            'err_Kprime_0':.2,
            'G_0': 109.0e9,
            'err_G_0':2e9,
            'Gprime_0': 1.2,
            'err_Gprime_0':0.1,
            'molar_mass': .4504,
            'n': 20,
            'Debye_0': 823.,
            'err_Debye_0':2.,
            'grueneisen_0': 1.05,
            'err_grueneisen_0':.06,
            'q_0': 1.9,
            'err_q_0':0.2,
            'eta_s_0': 2.4,
            'err_eta_s_0':0.1}


class almandine(Mineral):
    """
        garnet phase
        Stixrude & Lithgow-Bertelloni 2011 and references therein
        """
    def __init__(self,perturb=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 115.43e-6,
            'K_0': 174.0e9,
            'err_K_0':2e9,
            'Kprime_0': 4.9,
            'err_Kprime_0':.2,
            'G_0':96.0e9,
            'err_G_0':1e9,
            'Gprime_0': 1.4,
            'err_Gprime_0':0.1,
            'molar_mass': .4977,
            'n': 20,
            'Debye_0': 741.,
            'err_Debye_0':5.,
            'grueneisen_0': 1.06,
            'err_grueneisen_0':.06,
            'q_0': 1.4,
            'err_q_0':1.0,
            'eta_s_0': 2.1,
            'err_eta_s_0':1.0}

class mg_majorite(Mineral):
    """
        garnet phase
        Stixrude & Lithgow-Bertelloni 2011 and references therein
        """
    def __init__(self,perturb=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 114.32e-6,
            'K_0': 165.0e9,
            'err_K_0':3e9,
            'Kprime_0': 4.2,
            'err_Kprime_0':.3,
            'G_0': 85.0e9,#85.e9
            'err_G_0':2e9,
            'Gprime_0': 1.4, #1.4
            'err_Gprime_0':0.2,
            'molar_mass': .4016,
            'n': 20,
            'Debye_0': 822.,
            'err_Debye_0':4.,
            'grueneisen_0': 0.98,
            'err_grueneisen_0':.07,
            'q_0': 1.5,
            'err_q_0':0.5,
            'eta_s_0': 1.0,
            'err_eta_s_0':0.3}

class jd_majorite(Mineral):
    """
        garnet phase
        Stixrude & Lithgow-Bertelloni 2011 and references therein
        """
    def __init__(self,perturb=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 110.94e-6,
            'K_0': 177.0e9,
            'err_K_0':7e9,
            'Kprime_0': 4.1,
            'err_Kprime_0':1.0,
            'G_0': 125.0e9,
            'err_G_0':4e9,
            'Gprime_0': 1.4,
            'err_Gprime_0':0.5,
            'molar_mass': .4016,
            'n': 20,
            'Debye_0': 896.,
            'err_Debye_0':18.,
            'grueneisen_0': 1.01,
            'err_grueneisen_0':.3,
            'q_0': 1.4,
            'err_q_0':1.0,
            'eta_s_0': 3.3,
            'err_eta_s_0':1.}

class enstatite(Mineral):
    """
        orthopyroxene
        Stixrude & Lithgow-Bertelloni 2011 and references therein
        """
    def __init__(self,perturb=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 62.68e-6,
            'K_0': 107.0e9,
            'err_K_0':2e9,
            'Kprime_0': 7.0,
            'err_Kprime_0':.4,
            'G_0': 77.0e9,
            'err_G_0':1e9,
            'Gprime_0': 1.5,
            'err_Gprime_0':0.1,
            'molar_mass': .2008,
            'n': 10,
            'Debye_0': 812.,
            'err_Debye_0':4.,
            'grueneisen_0': 0.78,
            'err_grueneisen_0':.04,
            'q_0': 3.4,
            'err_q_0':0.4,
            'eta_s_0': 2.5,
            'err_eta_s_0':0.1}


class clinoenstatite(Mineral):
    """
        clinopyroxene
        Stixrude & Lithgow-Bertelloni 2011 and references therein
        """
    def __init__(self,perturb=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 62.50e-6,
            'K_0': 112.0e9,
            'err_K_0':10e9,
            'Kprime_0': 5.2,
            'err_Kprime_0':1.0,
            'G_0': 81.0e9,
            'err_G_0':10e9,
            'Gprime_0': 1.7,
            'err_Gprime_0':0.5,
            'molar_mass': .2008,
            'n': 10,
            'Debye_0': 805.,
            'err_Debye_0':10.,
            'grueneisen_0': 0.96,
            'err_grueneisen_0':0.3,
            'q_0': 1.5,
            'err_q_0':1.0,
            'eta_s_0': 1.7,
            'err_eta_s_0':0.5}

class hpclinoenstatite(Mineral):
    """
        HP-clinopyroxene
        Stixrude & Lithgow-Bertelloni 2011 and references therein
        """
    def __init__(self,perturb=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 60.76e-6,
            'K_0': 116.0e9,
            'err_K_0':1e9,
            'Kprime_0': 6.2,
            'err_Kprime_0':0.3,
            'G_0': 88.0e9,
            'err_G_0':1e9,
            'Gprime_0': 1.8,
            'err_Gprime_0':0.5,
            'molar_mass': .2008,
            'n': 10,
            'Debye_0': 824.,
            'err_Debye_0':7.,
            'grueneisen_0': 1.12,
            'err_grueneisen_0':0.05,
            'q_0': 0.2,
            'err_q_0':1.0,
            'eta_s_0': 2.1,
            'err_eta_s_0':0.5}

class diopside(Mineral):
    """
        clinopyroxene
        Stixrude & Lithgow-Bertelloni 2011 and references therein
        """
    def __init__(self,perturb=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 66.04e-6,
            'K_0': 112.0e9,
            'err_K_0':5e9,
            'Kprime_0': 5.2,
            'err_Kprime_0':1.8,
            'G_0': 67.0e9,
            'err_G_0':2e9,
            'Gprime_0': 1.4,
            'err_Gprime_0':0.5,
            'molar_mass': .2166,
            'n': 10,
            'Debye_0': 782.,
            'err_Debye_0':3.,
            'grueneisen_0': 0.96,
            'err_grueneisen_0':.05,
            'q_0': 1.5,
            'err_q_0':1.0,
            'eta_s_0': 1.6,
            'err_eta_s_0':1.0}



class ca_perovskite(Mineral):
    """
        Stixrude & Lithgow-Bertelloni 2011 and references therein
        """
    def __init__(self,perturb=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 27.45e-6,
            'K_0': 236.0e9,
            'err_K_0':4e9,
            'Kprime_0': 3.9,
            'err_Kprime_0':0.2,
            'G_0': 157.0e9,
            'err_G_0':12e9,
            'Gprime_0': 2.2,
            'err_Gprime_0':0.5,
            'molar_mass': .11616,
            'n': 5,
            'Debye_0': 796.,
            'err_Debye_0':44.,
            'grueneisen_0': 1.89,
            'err_grueneisen_0':.07,
            'q_0': 0.9,
            'err_q_0':1.6,
            'eta_s_0': 1.3,
            'err_eta_s_0':1.0}

class mgca_ferrite(Mineral):
    """
        Stixrude & Lithgow-Bertelloni 2011 and references therein
        """
    def __init__(self,perturb=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 36.18e-6,
            'K_0': 211.0e9,
            'err_K_0':1e9,
            'Kprime_0': 4.1,
            'err_Kprime_0':0.1,
            'G_0': 130.0e9,
            'err_G_0':1e9,
            'Gprime_0': 1.8,
            'err_Gprime_0':0.1,
            'molar_mass': .14226,
            'n': 7,
            'Debye_0': 838.,
            'err_Debye_0':16.,
            'grueneisen_0': 1.31,
            'err_grueneisen_0':.3,
            'q_0': 1.0,
            'err_q_0':1.0,
            'eta_s_0': 2.1,
            'err_eta_s_0':1.0}

class naca_ferrite(Mineral):
    """
        Stixrude & Lithgow-Bertelloni 2011 and references therein
        """
    def __init__(self,perturb=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 36.27e-6,
            'K_0': 158.0e9,
            'err_K_0':1e9,
            'Kprime_0': 4.3,
            'err_Kprime_0':0.1,
            'G_0': 121.0e9,
            'err_G_0':1e9,
            'Gprime_0': 2.1,
            'err_Gprime_0':0.1,
            'molar_mass': .14205,
            'n': 7,
            'Debye_0': 812.,
            'err_Debye_0':51.,
            'grueneisen_0': 1.17,
            'err_grueneisen_0':.3,
            'q_0': 1.0,
            'err_q_0':1.0,
            'eta_s_0': 1.6,
            'err_eta_s_0':1.0}

class corundum(Mineral):
    """
        akimotoite
        Stixrude & Lithgow-Bertelloni 2011 and references therein
        """
    def __init__(self,perturb=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 25.58e-6,
            'K_0': 253.0e9,
            'err_K_0':5e9,
            'Kprime_0': 4.3,
            'err_Kprime_0':0.2,
            'G_0': 163.0e9,
            'err_G_0':2e9,
            'Gprime_0': 1.6,
            'err_Gprime_0':0.1,
            'molar_mass': .10196,
            'n': 5,
            'Debye_0': 933.,
            'err_Debye_0':3.,
            'grueneisen_0': 1.32,
            'err_grueneisen_0':.4,
            'q_0': 1.3,
            'err_q_0':0.2,
            'eta_s_0': 2.8,
            'err_eta_s_0':0.2}

class mgakimotoite(Mineral):
    """
        akimotoite
        Stixrude & Lithgow-Bertelloni 2011 and references therein
        """
    def __init__(self,perturb=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 26.35e-6,
            'K_0': 211.0e9,
            'err_K_0':4e9,
            'Kprime_0': 5.6,
            'err_Kprime_0':0.8,
            'G_0': 132.0e9,
            'err_G_0':8e9,
            'Gprime_0': 1.6,
            'err_Gprime_0':0.5,
            'molar_mass': .10039,
            'n': 5,
            'Debye_0': 934.,
            'err_Debye_0':12.,
            'grueneisen_0': 1.19,
            'err_grueneisen_0':.13,
            'q_0': 2.3,
            'err_q_0':0.8,
            'eta_s_0': 2.8,
            'err_eta_s_0':1.0}
