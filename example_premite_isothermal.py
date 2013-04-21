import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
        sys.path.insert(1,os.path.abspath('..')) 

import burnman
from burnman import minerals

import pymc

seismic_model = burnman.seismic.prem() # pick from .prem() .slow() .fast() (see code/seismic.py)
number_of_points = 20 #set on how many depth slices the computations should be done
depths = np.linspace(750,2890, number_of_points)
seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)

print "preparations done"


def calc_velocities(ref_rho, ref_K, K_prime, ref_mu, mu_prime): 

    test = minerals.test_mineral()

    test.params['ref_V'] = 10.e-6
    test.params['molar_mass'] = ref_rho*test.params['ref_V']
    test.params['ref_K'] = ref_K
    test.params['K_prime'] = K_prime
    test.params['ref_mu'] = ref_mu
    test.params['mu_prime'] = mu_prime

    rock = burnman.composite( ( ( test, 1.0 ), ) )
    rock.set_method('bm')

    temperature = np.empty_like(seis_p)
    mat_rho, mat_vp, mat_vs, mat_vphi, mat_K, mat_mu = burnman.calculate_velocities(seis_p, temperature, rock)	

    return mat_rho, mat_vphi, mat_vs

def error(ref_rho, ref_K, K_prime, ref_mu, mu_prime):
    rho, vphi, vs = calc_velocities(ref_rho, ref_K, K_prime, ref_mu, mu_prime)

    vphi_chi = burnman.chi_factor(vphi, seis_vphi)
    vs_chi = burnman.chi_factor(vs, seis_vs)
    rho_chi = burnman.chi_factor(rho, seis_rho)

    return rho_chi+vphi_chi+vs_chi


# Priors on unknown parameters:
ref_rho = pymc.Uniform('ref_rho', lower=3300., upper=4500.)
ref_K = pymc.Uniform('ref_K', lower=200.e9, upper=300.e9)
K_prime = pymc.Uniform('K_prime', lower=3., upper=6.)
ref_mu = pymc.Uniform('ref_mu', lower=50.e9, upper=250.e9)
mu_prime = pymc.Uniform('mu_prime', lower=0., upper=3.)


minerr = 1e100

@pymc.deterministic
def theta(p1=ref_rho,p2=ref_K,p3=K_prime,p4=ref_mu,p5=mu_prime):
    global minerr
    if (p1<0 or p2<0 or p3<0 or p4<0 or p5 < 0):
        return 1e30
    try:
        e = error(p1,p2,p3,p4,p5)
        if (e<minerr):
            minerr=e
            print "best fit", e, "values:", p1,p2/1.e9,p3,p4/1.e9,p5
        return e
    except ValueError:
        return 1e20


sig = 1e-4
misfit = pymc.Normal('d',mu=theta,tau=1.0/(sig*sig),value=0,observed=True,trace=True)
model = [ref_rho, ref_K, K_prime, ref_mu, mu_prime, misfit]
things = ['ref_rho', 'ref_K', 'K_prime', 'ref_mu', 'mu_prime']

dbname = 'pickles/premite.pickle'
whattodo = ""

if len(sys.argv)<3:
    print "options:"
    print "run <dbname>"
    print "continue <dbname>"
    print "plot <dbname1> <dbname2> ..."
else:
    whattodo = sys.argv[1]
    dbname = sys.argv[2]
   
if whattodo=="run":
    S = pymc.MCMC(model, db='pickle', dbname=dbname)
    S.sample(iter=100, burn=0, thin=1)
    S.db.close()
    whattodo="continue"

if whattodo=="continue":
    n_runs = 100
    for l in range(0,n_runs):
        db = pymc.database.pickle.load(dbname)
        print "*** run=%d/%d, # samples: %d" % (l, n_runs, db.trace('ref_K').stats()['n'] )
        S = pymc.MCMC(model, db=db)
        S.sample(iter=500, burn=0, thin=1)
        S.db.close()
    whattodo="plot"

if whattodo=="plot":
    files=sys.argv[2:]
    print "files:",files

    b=10
    i=1

    for t in things:
        if t=='misfit':
            continue
	trace=[]
        print "trace:",t
        for filename in files:
            db = pymc.database.pickle.load(filename)
            newtrace=db.trace(t,chain=None).gettrace(burn=b,chain=None)
            if (trace!=[]):
                trace = np.append(trace, newtrace)
            else:
                trace=newtrace
                print "   adding ", newtrace.size, "burn = ",b
        print "  total size ", trace.size
        print "mean = ", trace.mean()
        for bin in [10,20,50,100]:
            hist,bin_edges=np.histogram(trace,bins=bin)
            a=np.argmax(hist)
            print "maxlike = ", bin_edges[a], bin_edges[a+1], (bin_edges[a]+bin_edges[a+1])/2.0

        pymc.Matplot.histogram(np.array(trace),t,rows=2,columns=(len(things)+1)/2,num=i)
        i=i+1

    plt.show()

if whattodo=="show":
	values = [float(i) for i in sys.argv[2:]]
	rho, vphi, vs = calc_velocities(values[0], values[1]*1.e9, values[2], values[3]*1.e9, values[4])

	plt.subplot(2,2,1)
	plt.plot(seis_p/1.e9,vs,color='r',linestyle='-',marker='^',markerfacecolor='r',markersize=4)
	plt.plot(seis_p/1.e9,seis_vs,color='k',linestyle='-',marker='v',markerfacecolor='k',markersize=4)
	plt.ylim([4, 8])
	plt.title("Vs (km/s)")

	plt.subplot(2,2,2)
	plt.plot(seis_p/1.e9,vphi,color='r',linestyle='-',marker='^',markerfacecolor='r',markersize=4)
	plt.plot(seis_p/1.e9,seis_vphi,color='k',linestyle='-',marker='v',markerfacecolor='k',markersize=4)
	plt.ylim([7, 12])
	plt.title("Vphi (km/s)")
    
    
	# plot density
	plt.subplot(2,2,3)
	plt.plot(seis_p/1.e9,rho,color='r',linestyle='-',marker='^',markerfacecolor='r',markersize=4,label='model 1')
	plt.plot(seis_p/1.e9,seis_rho,color='k',linestyle='-',marker='v',markerfacecolor='k',markersize=4,label='ref')
	plt.title("density (kg/m^3)")
	plt.legend(loc='upper left')
	plt.ylim([3, 7 ])
	plt.show()

print "done"


