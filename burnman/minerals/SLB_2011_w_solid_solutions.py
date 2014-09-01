# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""

SLB_2011
^^^^^^^^

Minerals from Stixrude & Lithgow-Bertelloni 2011 and references therein
Solid solutions from inv251010 of HeFESTo

"""

from burnman.processchemistry import ProcessChemistry
from burnman.mineral import Mineral
from burnman.solidsolution import SolidSolution

class c2c_pyroxene(SolidSolution):
    def __init__(self):
        # Name
        self.name='C2/c pyroxene'

        # Endmembers (C2/c is symmetric)
        base_material = [[hp_clinoenstatite(), '[Mg]2Si2O6'],[hp_clinoferrosilite(), '[Fe]2Si2O6']]

        # Interaction parameters (C2c is ideal)
        interaction_parameter=[]

        SolidSolution.__init__(self, base_material, interaction_parameter)



class ca_ferrite_structured_phase(SolidSolution):
    def __init__(self):
        # Name
        self.name='calcium ferrite structured phase'

        # Endmembers (CF is symmetric)
        base_material = [[mg_calcium_ferrite(), '[Mg]Al[Al]O4'],[fe_calcium_ferrite(), '[Fe]Al[Al]O4'],[na_calcium_ferrite(), '[Na]Al[Si]O4']]

        # Interaction parameters (CF is ideal)
        interaction_parameter=[]

        SolidSolution.__init__(self, base_material, interaction_parameter)

class ca_rich_perovskite(SolidSolution):
    def __init__(self):
        # Name
        self.name='calcium perovskite'

        # Endmembers (cpv is symmetric)
        base_material = [[ca_perovskite(), '[Ca][Ca]Si2O6'],[jadeitic_perovskite(), '[Na][Al]Si2O6']]

        # Interaction parameters (cpv is ideal)
        interaction_parameter=[]

        SolidSolution.__init__(self, base_material, interaction_parameter)

class clinopyroxene(SolidSolution):
    def __init__(self):
        # Name
        self.name='clinopyroxene'

        # Endmembers (cpx is symmetric)
        base_material = [[diopside(), '[Ca][Mg][Si]2O6'],[hedenbergite(), '[Ca][Fe][Si]2O6'],[clinoenstatite(), '[Mg][Mg][Si]2O6'],[ca_tschermaks_molecule(), '[Ca][Al][Si1/2Al1/2]2O6'],[jadeite(), '[Na][Al][Si]2O6']]

        # Interaction parameters
        interaction_parameter=[[0., 24.74e3, 26.e3, 24.3e3],[24.74e3, 0., 0.e3], [60.53136e3, 0.0], [10.e3]]

        SolidSolution.__init__(self, base_material, interaction_parameter)

class garnet(SolidSolution):
    def __init__(self):
        # Name
        self.name='garnet'

        # Endmembers (garnet is symmetric)
        base_material = [[pyrope(), '[Mg]3[Al][Al]Si3O12'],[almandine(), '[Fe]3[Al][Al]Si3O12'],[grossular(), '[Ca]3[Al][Al]Si3O12'],[mg_majorite(), '[Mg]3[Mg][Si]Si3O12'],[jadeitic_majorite(), '[Na2/3Al1/3]3[Al][Si]Si3O12']]
        # Interaction parameters
        interaction_parameter=[[0.0, 30.e3, 21.20278e3, 0.0],[0.0,0.0,0.0],[57.77596e3, 0.0],[0.0]]

        SolidSolution.__init__(self, base_material, interaction_parameter)


class akimotoite(SolidSolution):
    def __init__(self):
        # Name
        self.name='akimotoite/ilmenite'

        # Endmembers (ilmenite/akimotoite is symmetric)
        base_material = [[mg_akimotoite(), '[Mg][Si]O3'],[fe_akimotoite(), '[Fe][Si]O3'],[corundum(), '[Al][Al]O3']]
        # Interaction parameters
        interaction_parameter=[[0.0, 66.e3],[66.e3]]

        SolidSolution.__init__(self, base_material, interaction_parameter)

class ferropericlase(SolidSolution):
    def __init__(self):
        # Name
        self.name='magnesiowustite/ferropericlase'

        # Endmembers (ferropericlase is symmetric)
        base_material = [[periclase(), '[Mg]O'],[wuestite(), '[Fe]O']]
        # Interaction parameters
        interaction_parameter=[[13.e3]]

        SolidSolution.__init__(self, base_material, interaction_parameter)

class mg_fe_olivine(SolidSolution):
    def __init__(self):
        # Name
        self.name='olivine'

        # Endmembers (olivine is symmetric)
        base_material = [[forsterite(), '[Mg]2SiO4'],[fayalite(), '[Fe]2SiO4']]
        # Interaction parameters
        interaction_parameter=[[7.81322e3]]

        SolidSolution.__init__(self, base_material, interaction_parameter)

class orthopyroxene(SolidSolution):
    def __init__(self):
        # Name
        self.name='orthopyroxene'

        # Endmembers (orthopyroxene is symmetric)
        base_material = [[enstatite(), '[Mg][Mg][Si]SiO6'],[ferrosilite(), '[Fe][Fe][Si]SiO6'],[mg_tschermaks_molecule(), '[Mg][Al][Al]SiO6'],[orthodiopside(), '[Ca][Mg][Si]SiO6']]

        # Interaction parameters
        interaction_parameter=[[0.0, 0.0, 32.11352e3],[0.0, 0.0],[48.35316e3]]

        SolidSolution.__init__(self, base_material, interaction_parameter)

class plagioclase(SolidSolution):
    def __init__(self):
        # Name
        self.name='plagioclase'

        # Endmembers (plagioclase is symmetric)
        base_material = [[anorthite(), '[Ca][Al]2Si2O8'],[albite(), '[Na][Al1/2Si1/2]2Si2O8']]
        # Interaction parameters
        interaction_parameter=[[26.0e3]]

        SolidSolution.__init__(self, base_material, interaction_parameter)

class post_perovskite(SolidSolution):
    def __init__(self):
        # Name
        self.name='post-perovskite/bridgmanite'

        # Endmembers (post perovskite is symmetric)
        base_material = [[mg_post_perovskite(), '[Mg][Si]O3'],[fe_post_perovskite(), '[Fe][Si]O3'],[al_post_perovskite(), '[Al][Al]O3']]

        # Interaction parameters
        interaction_parameter=[[0.0, 60.0e3],[0.0]]

        SolidSolution.__init__(self, base_material, interaction_parameter)

class mg_fe_perovskite(SolidSolution):
    def __init__(self):
        # Name
        self.name='magnesium silicate perovskite/bridgmanite'

        # Endmembers (post perovskite is symmetric)
        base_material = [[mg_perovskite(), '[Mg][Si]O3'],[fe_perovskite(), '[Fe][Si]O3'],[al_perovskite(), '[Al][Al]O3']]

        # Interaction parameters
        interaction_parameter=[[0.0, 116.0e3],[0.0]]

        SolidSolution.__init__(self, base_material, interaction_parameter)

class mg_fe_ringwoodite(SolidSolution):
    def __init__(self):
        # Name
        self.name='ringwoodite'

        # Endmembers (post perovskite is symmetric)
        base_material = [[mg_ringwoodite(), '[Mg]2SiO4'],[fe_ringwoodite(), '[Fe]2SiO4']]

        # Interaction parameters
        interaction_parameter=[[9.34084e3]]

        SolidSolution.__init__(self, base_material, interaction_parameter)

class mg_fe_aluminous_spinel(SolidSolution):
    def __init__(self):
        # Name
        self.name='spinel-hercynite binary, fixed order'

        # Endmembers (post perovskite is symmetric)
        base_material = [[spinel(), '[Mg3/4Al1/4]4[Al7/8Mg1/8]8O16'],[hercynite(), '[Fe3/4Al1/4]4[Al7/8Fe1/8]8O16']]

        # Interaction parameters
        interaction_parameter=[[5.87646e3]]

        SolidSolution.__init__(self, base_material, interaction_parameter)

class mg_fe_wadsleyite(SolidSolution):
    def __init__(self):
        # Name
        self.name='wadsleyite'

        # Endmembers (post perovskite is symmetric)
        base_material = [[mg_wadsleyite(), '[Mg]2SiO4'],[fe_wadsleyite(), '[Fe]2SiO4']]

        # Interaction parameters
        interaction_parameter=[[16.74718e3]]

        SolidSolution.__init__(self, base_material, interaction_parameter)

class stishovite (Mineral):
    def __init__(self):
        self.params = {
            'equation_of_state': 'slb3',
            'V_0': 14.02e-6,
            'K_0': 314.0e9,
            'Kprime_0': 3.8,
            'G_0': 220.0e9,
            'Gprime_0': 1.9,
            'molar_mass': .0601,
            'n': 3,
            'Debye_0': 1108.,
            'grueneisen_0': 1.37,
            'q_0': 2.8,
            'eta_s_0': 4.6}

        self.uncertainties = {
             'err_K_0':8.e9,
             'err_Kprime_0':0.1,
             'err_G_0':12.e9,
             'err_Gprime_0':0.1,
             'err_Debye_0' : 13.,
             'err_grueneisen_0': .17,
             'err_q_0': 2.2,
             'err_eta_s_0' : 1.0
            }


class periclase (Mineral):
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
            'Debye_0': 767.,
            'grueneisen_0': 1.36,
            'q_0': 1.7, #1.7
            'eta_s_0': 2.8 } # 2.8

        self.uncertainties = {
        'err_K_0': 3.e9,
        'err_Kprime_0':.2,
        'err_G_0':1.0e9,
        'err_Gprime_0':.1,
        'err_Debye_0':9.,
        'err_grueneisen_0':.05,
        'err_q_0':.2,
        'err_eta_s_0':.2 }


class wuestite (Mineral):
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 12.26e-6,
            'K_0': 179.0e9,
            'Kprime_0': 4.9,
            'G_0': 59.0e9,
            'Gprime_0': 1.4,
            'molar_mass': .0718,
            'n': 2,
            'Debye_0': 454.,
            'grueneisen_0': 1.53,
            'q_0': 1.7, #1.7
            'eta_s_0': -0.1 } #

        self.uncertainties = {
            'err_K_0':1.e9,
            'err_Kprime_0':.2,
            'err_G_0':1.e9,
            'err_Gprime_0':.1,
            'err_Debye_0':21.,
            'err_grueneisen_0':.13,
            'err_q_0':1.0,
            'err_eta_s_0':1.0}

class mg_perovskite(Mineral):
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 24.45e-6,
            'K_0': 251.0e9,
            'Kprime_0': 4.1,
            'G_0': 173.0e9,
            'Gprime_0': 1.7,
            'molar_mass': .1000,
            'n': 5,
            'Debye_0': 905.,
            'grueneisen_0': 1.57,
            'q_0': 1.1,
            'eta_s_0': 2.6 } #2.6

        self.uncertainties = {
            'err_K_0': 3.e9,
            'err_Kprime_0': 0.1,
            'err_G_0': 2.e9,
            'err_Gprime_0' : 0.0,
            'err_Debye_0': 5.,
            'err_grueneisen_0':.05,
            'err_q_0': .3,
            'err_eta_s_0':.3}


class fe_perovskite(Mineral):
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 25.49e-6,
            'K_0': 272.0e9,
            'Kprime_0': 4.1,
            'G_0': 133.0e9,
            'Gprime_0': 1.4,
            'molar_mass': .1319,
            'n': 5,
            'Debye_0': 871.,
            'grueneisen_0': 1.57,
            'q_0': 1.1,
            'eta_s_0': 2.3 } #2.3

        self.uncertainties = {
            'err_K_0':40e9,
            'err_Kprime_0':1.,
            'err_G_0':40e9,
            'err_Gprime_0':0.0,
            'err_Debye_0':26.,
            'err_grueneisen_0':.3,
            'err_q_0':1.0,
            'err_eta_s_0':1.0}

# Mineral aliases
mgc2 = hp_clinoenstatite
fec2 = hp_clinoferrosilite
mgcf = mg_calcium_ferrite
fecf = fe_calcium_ferrite
nacf = na_calcium_ferrite
capv = calcium_perovskite
jdpv = jadeitic_perovskite
di = diopside
he = hedenbergite
cen = clinoenstatite
cats = ca_tschermaks_molecule
jd = jadeite
py = pyrope
al = almandine
gr = grossular
mgmj = mg_majorite
jdmj = jadeitic_majorite
mgil = mg_akimotoite
feil = fe_akimotoite
co = corundum
pe = periclase
wu = wuestite
fo = forsterite
fa = fayalite
en = enstatite
fs = ferrosilite
mgts = mg_tschermaks_molecule
odi = orthodiopside
an = anorthite
ab = albite
mppv = mg_post_perovskite
fppv = fe_post_perovskite
appv = al_post_perovskite
mgpv = mg_perovskite
mg_bridgmanite = mg_perovskite
fepv = fe_perovskite
fe_bridgmanite = fe_perovskite
alpv = al_perovskite
mgri = mg_ringwoodite
feri = fe_ringwoodite
sp = spinel
hc = hercynite
mgwa = mg_wadsleyite
fewa = fe_wadsleyite


# Solid solution aliases
c2c = c2c_pyroxene
cf = calcium_ferrite_structured_phase
cpv = ca_rich_perovskite
cpx = clinopyroxene
gt = garnet
il = akimotoite
ilmenite_group = akimotoite
mw = ferropericlase
magnesiowuestite = ferropericlase
ol = mg_fe_olivine
opx = orthopyroxene
plag = plagioclase
ppv = post_perovskite
pv = mg_fe_perovskite
mg_fe_bridgmanite = mg_fe_perovskite
mg_fe_silicate_perovskite = mg_fe_perovskite
ri = mg_fe_ringwoodite
spinel_group=mg_fe_aluminous_spinel
wa = mg_fe_wadsleyite
spinelloid_III = mg_fe_wadsleyite
