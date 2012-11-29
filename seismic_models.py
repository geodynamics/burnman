import numpy
import bisect
import matplotlib.pyplot as pyplot
from tools import *
import seismo_in as seis

# this loads different seismic models
# available models: 'prem', 'ref_fast', 'ref_slow'

# pressure: in GPa
# return: V_p, V_s in km/s
    #def prem_V(pressure):
#return lookup_(pressure, 3), lookup_(pressure, 4)

table = []

def prem_attenuation(pressure):
    return lookup_('prem',pressure, 5), lookup_('prem',pressure, 6)

def prem_depth(pressure):
    return lookup_('prem',pressure, 0)

def prem_density(pressure):
    return lookup_('prem',pressure, 2)

def lookup_depthbounds(name):
    get_model(name)
    return min(table[:,0]), max(table[:,0])

#return linear interpolated column colidx
def lookup_(name,pressure, colidx):
    get_model(name)
    idx = bisect.bisect_left(table[:,1], pressure) - 1
    if (idx < 0):
        return table[0][colidx]
    elif (idx < len(table[:,1])-1):
        return linear_interpol(pressure, table[idx,1], table[idx+1,1], table[idx][colidx], table[idx+1][colidx])
    else:
        return table[idx][colidx]

#return pressure for a given radius from a lookup table
def lookup_pressure(name,radius):
    get_model(name)
    idx = (bisect.bisect_left(table[:,0], radius)) - 1    
    if (idx < 0):
        return table[1,1]
    elif (idx < len(table[:,0])-1):
        return linear_interpol(radius, table[idx,0], table[idx+1,0], table[idx,1], table[idx+1,1])
    else:
        return table[idx,1]

def get_model(name):
    # Models implemented so far: 'prem', 'ref_fast', 'ref_slow'
    # read model in global variable table. Data is structure as:
    #radius pressure density V_p V_s Q_K Q_mu
    global table
    table = read_table("data/prem_withQ.txt")
    table = sort_table(table, 1)
    table = numpy.array(table)
    table[:,0]=6371.0e3-table[:,0] # converts radius to depth
    table_p = table[:,1]
    if (name != "prem"):
        #Read in fast and slow reference models
        if (name=="ref_fast"):
            table_news = read_table("data/swave_fast.txt")
            table_newp = read_table("data/pwave_fast.txt")
        elif (name=="ref_slow"):       
            table_news = read_table("data/swave_slow.txt")
            table_newp = read_table("data/pwave_slow.txt")
        table_news = sort_table(table_news,1)
        table_news = numpy.array(table_news)
        table_newp = sort_table(table_newp,1)
        table_newp = numpy.array(table_newp)   
        idx_min = bisect.bisect_left(table[:,0], min(table_news[:,0])) - 1
        idx_max = bisect.bisect_left(table[:,0], max(table_news[:,0])) 
        table = cut_table(table,idx_min,idx_max)
        table = sort_table(table, 1)
        table = numpy.array(table)
        table[:,3]=table_newp[:,1]
        table[:,4]=table_news[:,1]
    return

                         
# run seismic_models.py to plot the different seismic models available
if __name__ == "__main__":
    get_model('prem')
    pyplot.subplot(2,2,1)
    p1,=pyplot.plot(table[:,1],table[:,3]/1.0e3,'-k',linewidth=2.0)
    pyplot.subplot(2,2,2)
    pyplot.plot(table[:,1],table[:,4]/1.0e3,'-k',linewidth=2.0)
    get_model('ref_fast')
    pyplot.subplot(2,2,1)
    p2,=pyplot.plot(table[:,1],table[:,3]/1.0e3,'-b',linewidth=2.0)
    pyplot.subplot(2,2,2)
    pyplot.plot(table[:,1],table[:,4]/1.0e3,'-b',linewidth=2.0)  
    get_model('ref_slow')
    pyplot.subplot(2,2,1)
    p3,=pyplot.plot(table[:,1],table[:,3]/1.0e3,'-r',linewidth=2.0)
    pyplot.subplot(2,2,2)
    pyplot.plot(table[:,1],table[:,4]/1.0e3,'-r',linewidth=2.0)
    
    pyplot.subplot(2,2,1)
    pyplot.legend([p1,p2,p3],["PREM","ref_fast","ref_slow"], loc=4)
    pyplot.xlabel('pressure (Gpa)')
    pyplot.ylabel('velocity (km/s)')
    pyplot.title('P wave velocity')
    pyplot.subplot(2,2,2)
    pyplot.xlabel('pressure (Gpa)')
    pyplot.ylabel('velocity (km/s)')
    pyplot.title('S wave velocity')
    pyplot.show()
