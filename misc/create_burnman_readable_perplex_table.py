from subprocess import Popen, PIPE, STDOUT

'''
# Properties
    1 - Specific Enthalpy (J/m3)                                    
X    2 - Density (kg/m3)                                             
    3 - Specific heat capacity (J/K/m3)                             
X    4 - Expansivity (1/K, for volume)                               
X    5 - Compressibility (1/bar, for volume)                         
    6 - Composition (Mol or Wt%) of the system                      
    7 - Mode (Vol, Mol, or Wt proportion) of a phase                
    8 - Composition (Mol or Wt%) of a solution phase                
    9 - Grueneisen thermal ratio                                    
X   10 - Adiabatic bulk modulus (bar)                                
X   11 - Adiabatic shear modulus (bar)                               
X   12 - Sound velocity (km/s)                                       
X   13 - P-wave velocity (Vp, km/s)                                  
X   14 - S-wave velocity (Vs, km/s)                                  
   15 - Vp/Vs                                                       
   16 - Specific entropy (J/K/m3)                                   
X   17 - Entropy (J/K/kg)                                            
X   18 - Enthalpy (J/kg)                                             
X   19 - Heat Capacity (J/K/kg)                                      
   20 - Specific mass of a phase (kg/m3-system)                     
   21 - Poisson ratio                                               
X   22 - Molar Volume (J/bar)                                        
   23 - Dependent potentials (J/mol, bar, K)                        
   24 - Assemblage Index                                            
   25 - Modes of all phases                                         
   26 - Sound velocity T derivative (km/s/K)                        
   27 - P-wave velocity T derivative (km/s/K)                       
   28 - S-wave velocity T derivative (km/s/K)                       
   29 - Adiabatic bulk modulus T derivative (bar/K)                 
   30 - Shear modulus T derivative (bar/K)                          
   31 - Sound velocity P derivative (km/s/bar)                      
   32 - P-wave velocity P derivative (km/s/bar)                     
   33 - S-wave velocity P derivative (km/s/bar)                     
   34 - Adiabatic bulk modulus P derivative (unitless)              
   35 - Shear modulus P derivative (unitless)                       
   36 - All phase &/or system properties                            
   37 - Absolute amount (Vol, Mol, or Wt) of a phase                
   38 - Multiple property output                                    
   39 - Heat capacity ratio (Cp/Cv)
'''

werami_path = '../programs_6.6.8_20140304/werami'
project = 'in23'
n_pressures = 11
n_temperatures = 11
change_PT_range = 'y'
P_range = [10.e9, 110.e9] # Pa
T_range = [300., 3000.] # K

if change_PT_range == 'n':
    str2 = 'n\n'
else:
    str2 = 'y\n{0} {1}\n{2} {3}\n'.format(P_range[0]/1.e5, P_range[1]/1.e5,
                                          T_range[0], T_range[1])
    
stdin='{0:s}\n' \
    '2\n' \
    '2\n' \
    'n\n' \
    '4\n' \
    'n\n' \
    '5\n' \
    'n\n' \
    '10\n' \
    'n\n' \
    '11\n' \
    'n\n' \
    '12\n' \
    'n\n' \
    '13\n' \
    'n\n' \
    '14\n' \
    'n\n' \
    '17\n' \
    'n\n' \
    '18\n' \
    'n\n' \
    '19\n' \
    'n\n' \
    '22\n' \
    'n\n' \
    '0\n' \
    '{1:s}' \
    '{2:d} {3:d}\n' \
    '0'.format(project, str2, n_pressures, n_temperatures)


p = Popen([werami_path], stdout=PIPE, stdin=PIPE, stderr=STDOUT)

print('Working on creating {0}x{1} P-T table file using werami. Please wait.'.format(n_pressures, n_temperatures))
stdout = p.communicate(input=stdin)[0]
print(stdout)
print('Processing complete')
