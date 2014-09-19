import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..'))


import os, sys
from time import time
import pymc
import math
import cProfile
from scipy.stats import norm
import matplotlib.mlab as mlab


if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..'))
import burnman
from burnman import minerals

seismic_model = burnman.seismic.PREM() # pick from .prem() .slow() .fast() (see code/seismic.py)
number_of_points = 10 #set on how many depth slices the computations should be done
depths = np.linspace(1000e3,2500e3, number_of_points)
seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)


temperature = burnman.geotherm.brown_shankland(seis_p)


print "preparations done"


# rewrite as function of f_pv and fe_content
def calc_velocities(amount_perovskite, fe_pv, fe_pc):
    print amount_perovskite,fe_pv,fe_pc
    method = 'slb3' #slb3|slb2|mgd3|mgd2
    pv=minerals.SLB_2011.mg_fe_perovskite(fe_pv)
    pc=minerals.SLB_2011.ferropericlase(fe_pc)
    rock = burnman.Composite( [amount_perovskite, 1.0-amount_perovskite],[ pv,pc ] )
    
    
    rock.set_method(method)
    
    mat_rho, mat_vp, mat_vs, mat_vphi, mat_K, mat_G = burnman.velocities_from_rock(rock,seis_p, temperature)
    return mat_vp, mat_vs, mat_rho

def error(amount_perovskite, fe_pv, fe_pc):
	mat_vp, mat_vs, mat_rho = calc_velocities(amount_perovskite, fe_pv, fe_pc)
	
	vs_err = burnman.l2(depths, mat_vs, seis_vs)
	vp_err = burnman.l2(depths, mat_vp, seis_vp)
	den_err = burnman.l2(depths, mat_rho, seis_rho)
    
	return vs_err + vp_err + den_err


# Priors on unknown parameters:
# make flat distibutions
fraction_pv = pymc.Uniform('fraction_pv', 0.0, 1.0)
fe_pv = pymc.Uniform('fe_pv', 0.0, 0.5)
fe_pc= pymc.Uniform('fe_pc', 0.0, 0.5)



minerr = 1e100

#(plot=False)
@pymc.deterministic
def theta(p1=fraction_pv,p2=fe_pv,p3=fe_pc):
	global minerr
    
	if (p1<0. or p2<0. or p3<0. or p1>1.0 or p2>1.0 or p3>1.0):
         return 1e30
    
	try:
		e = error(p1,p2,p3)
		if (e<minerr):
			minerr=e
			print "best fit", e, "values:", p1,p2,p3
		return e
	except ValueError:
		return 1e20#float("inf")


sig = 1e-4
misfit = pymc.Normal('d',mu=theta,tau=1.0/(sig*sig),value=0,observed=True,trace=True)
model = [fraction_pv, fe_pv, fe_pc, misfit]
things = ['fraction_pv', 'fe_pv','fe_pc','misfit']

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
    n_runs = 200
    for l in range(0,n_runs):
        db = pymc.database.pickle.load(dbname)
        print "*** run=%d/%d, # samples: %d" % (l, n_runs, db.trace('fraction_pv').stats()['n'] )
        S = pymc.MCMC(model, db=db)
        S.sample(iter=1000, burn=100, thin=1)
        S.db.close()

if whattodo=="plot":
	files=sys.argv[2:]
	print "files:",files
        
        b=0#10000 # burn number
        i=1

        for t in things:
            if t=='misfit':
                continue
            trace=[]
            print "trace:",t
            for filename in files:
                db = pymc.database.pickle.load(filename)
                newtrace=db.trace(t,chain=None).gettrace(burn=b,chain=None)
                print newtrace.size
                if (trace!=[]):
                    trace = np.append(trace, newtrace)
                else:
                    trace=newtrace
                print "   adding ", newtrace.size, "burn = ",b
            print "  total size ", trace.size
            print "mean = ", trace.mean()
            print trace[0:100]
            for bin in [10,20,50,100]:
                hist,bin_edges=np.histogram(trace,bins=bin)
                a=np.argmax(hist)
                print "maxlike = ", bin_edges[a], bin_edges[a+1], (bin_edges[a]+bin_edges[a+1])/2.0
            
            
            (mu, sigma) = norm.fit(np.array(trace))
            print "mu, sigma: %e %e" % (mu, sigma)
            plt.subplot(2,len(things)/2,i)
            n, bins, patches = plt.hist(np.array(trace), 60, normed=1, facecolor='green', alpha=0.75)
            y = mlab.normpdf( bins, mu, sigma)
            l = plt.plot(bins, y, 'r--', linewidth=2)
            plt.title("%s, mean: %.3e, std dev.: %.3e" % (t,mu,sigma),fontsize='small')
            
            #pymc.Matplot.histogram(np.array(trace),t,rows=2,columns=len(things)/2,num=i)
            
            
            
            i=i+1
            
        plt.savefig("example_inv_big_pv.png")
        plt.show()

if whattodo=="test":
    db = pymc.database.pickle.load(dbname)
    S = pymc.MCMC(model, db=db)
    
    for t in things:
        print db.trace(t).stats()
    
    print "means:"
    for t in things:
	    print t,"\t",db.trace(t).stats()['mean']
    
    print "#samples: %d" % db.trace('mg_pv_K').stats()['n']
    
    pymc.raftery_lewis(S, q=0.025, r=0.01)
    
    b = 1
    t = 1
    
    scores = pymc.geweke(S, intervals=20)
    
    pymc.Matplot.trace(db.trace('deviance',chain=None).gettrace(burn=1000,thin=t,chain=None),'deviance',rows=2,columns=9,num=1)
    
    
    pymc.Matplot.trace(db.trace('fraction_pv',chain=None).gettrace(thin=t,chain=None),'fraction_pv',rows=2,columns=9,num=2)
    pymc.Matplot.histogram(np.array(db.trace('mg_pv_K',chain=None).gettrace(burn=b,chain=None)),'fraction_pv',rows=2,columns=9,num=11)
    
    pymc.Matplot.trace(db.trace('fe_pv',chain=None).gettrace(thin=t,chain=None),'fe_pv',rows=2,columns=9,num=3)
    pymc.Matplot.histogram(np.array(db.trace('fe_pv',chain=None).gettrace(burn=b,chain=None)),'fe_pv',rows=2,columns=9,num=12)
    
    pymc.Matplot.trace(db.trace('fe_pc',chain=None).gettrace(thin=t,chain=None),'fe_pv',rows=2,columns=9,num=4)
    pymc.Matplot.histogram(np.array(db.trace('fe_pv',chain=None).gettrace(burn=b,chain=None)),'fe_pv',rows=2,columns=9,num=13)
    
    
    
    plt.show()


if whattodo=="show":
	values = [float(i) for i in sys.argv[2:]]
	mat_vp, mat_vs, mat_rho = calc_velocities(values[0], values[1], values[2])
    
	plt.subplot(2,2,1)
	plt.plot(seis_p/1.e9,mat_vs/1.e3,color='r',linestyle='-',marker='^',markerfacecolor='r',markersize=4)
	plt.plot(seis_p/1.e9,seis_vs/1.e3,color='k',linestyle='-',marker='v',markerfacecolor='k',markersize=4)
	plt.ylim([4, 8])
	plt.title("Vs (km/s)")
    
    
	# plot Vphi
	plt.subplot(2,2,2)
	plt.plot(seis_p/1.e9,mat_vp/1.e3,color='r',linestyle='-',marker='^',markerfacecolor='r',markersize=4)
	plt.plot(seis_p/1.e9,seis_vp/1.e3,color='k',linestyle='-',marker='v',markerfacecolor='k',markersize=4)
	plt.ylim([10, 14])
	plt.title("Vp (km/s)")
    
	# plot density
	plt.subplot(2,2,3)
	plt.plot(seis_p/1.e9,mat_rho/1.e3,color='r',linestyle='-',marker='^',markerfacecolor='r',markersize=4,label='model 1')
	plt.plot(seis_p/1.e9,seis_rho/1.e3,color='k',linestyle='-',marker='v',markerfacecolor='k',markersize=4,label='ref')
	plt.title("density (kg/m^3)")
	plt.legend(loc='upper left')
	plt.ylim([4, 8])
	#plt.savefig("output_figures/example_inv_big_pv_show.png")
	plt.show()


if whattodo=="profile2":
	#run with:
	#python -m cProfile -o output.pstats example_inv_big_pv.py profile2 1
	#gprof2dot.py -f pstats output.pstats | dot -Tpng -o output.png
	[error(235.654790318e9, 3.87724833477, 165.45907725e9, 1.61618366689, 273.888690109e9, 4.38543140228, 306.635371217e9, 1.42610871557) for i in range(0,10)]

if whattodo=="profile":
	#just run normally
	cProfile.run("[error(235.654790318e9, 3.87724833477, 165.45907725e9, 1.61618366689, 273.888690109e9, 4.38543140228, 306.635371217e9, 1.42610871557) for i in range(0,10)]")



