from subprocess import Popen, PIPE, STDOUT
import argparse

'''  
This python file uses PerpleX's werami software 
to output a table file containing the following
material properties.            
    2 - Density (kg/m3)                           
    4 - Expansivity (1/K, for volume)                               
    5 - Compressibility (1/bar, for volume)                     
   10 - Adiabatic bulk modulus (bar)                                
   11 - Adiabatic shear modulus (bar)                               
   12 - Sound velocity (km/s)                                       
   13 - P-wave velocity (Vp, km/s)                                  
   14 - S-wave velocity (Vs, km/s)                               
   17 - Entropy (J/K/kg)                                            
   18 - Enthalpy (J/kg)                                             
   19 - Heat Capacity (J/K/kg)                              
   22 - Molar Volume (J/bar)        
'''

parser = argparse.ArgumentParser(description='Call werami to create a burnman-readable tab file.')

parser.add_argument('--werami_path', metavar='path', type=str, nargs='+', required=True,
                    help='The path to werami')
parser.add_argument('--project', metavar='project name', type=str, nargs='+', required=True,
                    help='The name of the project file (without the suffix)')
parser.add_argument('--n_pressures', type=int, nargs=1, required=True,
                    help='The number of pressure steps in the grid')
parser.add_argument('--n_temperatures', type=int, nargs=1, required=True,
                    help='The number of pressure steps in the grid')
parser.add_argument('--P_range', type=float, nargs=2,
                    help='Minimum and maximum values of pressure (Pa; optional)')
parser.add_argument('--T_range', type=float, nargs=2,
                    help='Minimum and maximum values of temperature (K; optional)')

args = parser.parse_args()

print('Working on creating {0}x{1} P-T table file using werami. Please wait.\n'.format(args.n_pressures[0], args.n_temperatures[0]))

try:
    str2 = 'y\n{0} {1}\n{2} {3}\n'.format(args.P_range[0]/1.e5, args.P_range[1]/1.e5,
                                          args.T_range[0], args.T_range[1])
except:
    print('Keeping P-T range the same as the original project range.\n'
          'If you wish to change the range, you will need to add the following command line arguments:\n'
          '--P_range [P_min] [P_max] (in Pa)\n'
          '--T_range [T_min] [T_max] (in K)\n')
    str2 = 'n\n'
    
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
    '0'.format(args.project[0], str2, args.n_pressures[0], args.n_temperatures[0])


p = Popen(args.werami_path, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
stdout = p.communicate(input=stdin)[0]
print(stdout)
print('Processing complete')
