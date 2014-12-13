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



number_of_points = 5 #set on how many depth slices the computations should be done



# velocity constraints from seismology
seismic_model = burnman.seismic.PREM() # pick from .prem() .slow() .fast() (see code/seismic.py)
depths = np.linspace(1000e3,2500e3, number_of_points)
seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)


seismic_model2 = burnman.seismic.PREM() # pick from .prem() .slow() .fast() (see code/seismic.py)
depths = np.linspace(1000e3,2500e3, number_of_points)
seis_p2, seis_rho2, seis_vp2, seis_vs2, seis_vphi2 = seismic_model2.evaluate_all_at(depths)

'''
print seis_vs-seis_vs2
plt.plot(depths,seis_vs,'r')
plt.plot(depths,seis_vs2,'k')
plt.plot(depths,seis_vp,'r')
plt.plot(depths,seis_vp2,'k')
#plt.plot(depths,seis_rho)
#plt.plot(depths,seis_rho2)
plt.show()
stophere
'''

seis_G= seis_vs**2.*seis_rho
seis_K= seis_vphi**2*seis_rho

'''

#velocities from known material (still named seis*)
amount_perovskite = 0.8
pv=minerals.SLB_2011.mg_fe_perovskite(.06)
pc=minerals.SLB_2011.ferropericlase(.2)
rock = burnman.Composite( [amount_perovskite, 1.0-amount_perovskite],[ pv,pc ] )
pressures = np.linspace(25e9,130e9,number_of_points)
temperature = burnman.geotherm.brown_shankland(pressures)
rock.set_method('slb3')
print "Calculations are done for:"
rock.debug_print()
seis_p=pressures
seis_rho, seis_vp, seis_vs, seis_vphi, seis_K, seis_G = \
    burnman.velocities_from_rock(rock, pressures, temperature, \
                                 burnman.averaging_schemes.VoigtReussHill())
depths= burnman.depths_for_rock(rock, pressures, temperature, \
                                 burnman.averaging_schemes.VoigtReussHill())

'''

print "preparations done"




# Priors on unknown parameters:
fraction_pv = pymc.Uniform('fraction_pv', 0.0, 1.0)
fe_pv = pymc.Uniform('fe_pv', 0.0, 0.5)
fe_pc= pymc.Uniform('fe_pc', 0.0, 0.5)
T0 = pymc.Normal('T0',mu=1500.,tau=1./500.)
def calc_all_velocities(fraction_pv, fe_pv, fe_pc,T0):
        method = 'slb3' #slb3|slb2|mgd3|mgd2
        pv=minerals.other.mg_fe_perovskite(fe_pv)
        pc=minerals.other.ferropericlase(fe_pc)
        rock = burnman.Composite( [fraction_pv, 1.0-fraction_pv],[ pv,pc ] )
        temperature = burnman.geotherm.adiabatic(seis_p,T0, rock)


        mat_rho, mat_vp, mat_vs, mat_vphi, mat_K, mat_G = burnman.velocities_from_rock(rock,seis_p, temperature)
        return mat_vp, mat_vs, mat_rho, mat_vphi, mat_K, mat_G

def nrmse(funca,funcb):
    """
    Normalized root mean square error for one profile
    :type funca: list of arrays of float
    :param funca: array calculated values
    :type funcb: list of arrays of float
    :param funcb: array of values (observed or calculated) to compare to

    :returns: RMS error
    :rtype: array of floats

    """
    diff=np.array(funca-funcb)
    diff=diff*diff
    rmse=np.sqrt(np.sum(diff)/len(diff))
    nrmse=rmse/(np.max(funcb)-np.min(funcb))
    return nrmse

def error(fraction_pv, fe_pv, fe_pc,T0):
    if True:
        if fraction_pv>0. and fraction_pv<1.0 and fe_pv>0. and fe_pv<1.0 and fe_pc>0. and fe_pc<1.0 and T0>300.:
            mat_vp, mat_vs, mat_rho, mat_vphi, mat_K, mat_G = calc_all_velocities(fraction_pv, fe_pv, fe_pc,T0)

            misfit = nrmse(mat_rho,seis_rho)+nrmse(mat_K,seis_K)+nrmse(mat_G,seis_G)
            return misfit
        else:
            return 1.e30

    #except:
    #    return 1e30

# rewrite as function of f_pv and fe_content
@pymc.deterministic
def calc_misfit(fraction_pv=fraction_pv, fe_pv=fe_pv, fe_pc=fe_pc,T0=T0):
    return error(fraction_pv,fe_pv,fe_pc,T0)


#sig = 1e-2
#misfit = pymc.Normal('d',mu=theta,tau=1.0/(sig*sig),value=0.,observed=True,trace=True)
#sig = pymc.Uniform('sig', 0.0, 100.0, value=1.)
sig=1.e-2 # Some sorts of error
obs = pymc.Normal('d',mu=calc_misfit,tau=1.0/(sig*sig),value=0,observed=True,trace=True)
model = [fraction_pv, fe_pv, fe_pc, T0,obs]
things = ['fraction_pv', 'fe_pv','fe_pc','T0']

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
    #pymc.MAP(model).fit() # Find minimum to start search from
    S = pymc.MCMC(model, db='pickle', dbname=dbname)
    S.sample(iter=100, burn=0, thin=1)
    S.db.close()
    whattodo="continue"

if whattodo=="continue":
    n_runs = 100
    for l in range(0,n_runs):
        db = pymc.database.pickle.load(dbname)
        print "*** run=%d/%d, # samples: %d" % (l, n_runs, db.trace('fraction_pv').stats()['n'] )
        S = pymc.MCMC(model, db=db)
        #S.sample(iter=100, burn=10, thin=1)
        S.sample(iter=10, burn=0, thin=1) # Search space for 100000 acceptable steps, forget first 1000 and save every 10.
        S.db.close()

if whattodo=="plot":
	files=sys.argv[2:]
	print "files:",files

        b=0#10000 # burn number
        i=1

        for t in things:
            print t
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
            print trace[0:100]
            for bin in [10,20,50,100]:
                hist,bin_edges=np.histogram(trace,bins=bin)
                a=np.argmax(hist)
                print "maxlike = ", bin_edges[a], bin_edges[a+1], (bin_edges[a]+bin_edges[a+1])/2.0


            (mu, sigma) = norm.fit(np.array(trace))
            print "mu, sigma: %e %e" % (mu, sigma)
            plt.subplot(2,(len(things)+1)/2,i)
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

    print "#samples: %d" % db.trace('fraction_pv').stats()['n']

    pymc.raftery_lewis(S, q=0.025, r=0.01)

    b = 1
    t = 1

    scores = pymc.geweke(S, intervals=20)
    print scores
    pymc.Matplot.trace(db.trace('deviance',chain=None).gettrace(burn=1000,thin=t,chain=None),'deviance',rows=2,columns=4,num=1)


    pymc.Matplot.trace(db.trace('fraction_pv',chain=None).gettrace(thin=t,chain=None),'fraction_pv',rows=2,columns=4,num=2)
    pymc.Matplot.histogram(np.array(db.trace('fraction_pv',chain=None).gettrace(burn=b,chain=None)),'fraction_pv',rows=2,columns=4,num=6)

    pymc.Matplot.trace(db.trace('fe_pv',chain=None).gettrace(thin=t,chain=None),'fe_pv',rows=2,columns=4,num=3)
    pymc.Matplot.histogram(np.array(db.trace('fe_pv',chain=None).gettrace(burn=b,chain=None)),'fe_pv',rows=2,columns=4,num=7)

    pymc.Matplot.trace(db.trace('fe_pc',chain=None).gettrace(thin=t,chain=None),'fe_pc',rows=2,columns=4,num=4)
    pymc.Matplot.histogram(np.array(db.trace('fe_pc',chain=None).gettrace(burn=b,chain=None)),'fe_pc',rows=2,columns=4,num=8)



    plt.show()


if whattodo=="show":
	values = [float(i) for i in sys.argv[2:]]
	#mat_vp, mat_vs, mat_rho = calc_velocities(values[0], values[1], values[2])
        mat_vp,mat_vs, mat_rho, mat_vphi=calc_all_velocities(values[0], values[1], values[2])
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
	[error(0.5,0.1,0.2) for i in range(0,1000)]

if whattodo=="profile":
	#just run normally
        print error(0.5,0.1,0.2)
	cProfile.run("[error(0.5,0.1,0.2) for i in range(0,1000)]")import os, sys, numpy as np, matplotlib.pyplot as plt
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



number_of_points = 5 #set on how many depth slices the computations should be done



# velocity constraints from seismology
seismic_model = burnman.seismic.PREM() # pick from .prem() .slow() .fast() (see code/seismic.py)
depths = np.linspace(1000e3,2500e3, number_of_points)
seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)


seismic_model2 = burnman.seismic.PREM() # pick from .prem() .slow() .fast() (see code/seismic.py)
depths = np.linspace(1000e3,2500e3, number_of_points)
seis_p2, seis_rho2, seis_vp2, seis_vs2, seis_vphi2 = seismic_model2.evaluate_all_at(depths)

'''
print seis_vs-seis_vs2
plt.plot(depths,seis_vs,'r')
plt.plot(depths,seis_vs2,'k')
plt.plot(depths,seis_vp,'r')
plt.plot(depths,seis_vp2,'k')
#plt.plot(depths,seis_rho)
#plt.plot(depths,seis_rho2)
plt.show()
stophere
'''

seis_G= seis_vs**2.*seis_rho
seis_K= seis_vphi**2*seis_rho

'''

#velocities from known material (still named seis*)
amount_perovskite = 0.8
pv=minerals.SLB_2011.mg_fe_perovskite(.06)
pc=minerals.SLB_2011.ferropericlase(.2)
rock = burnman.Composite( [amount_perovskite, 1.0-amount_perovskite],[ pv,pc ] )
pressures = np.linspace(25e9,130e9,number_of_points)
temperature = burnman.geotherm.brown_shankland(pressures)
rock.set_method('slb3')
print "Calculations are done for:"
rock.debug_print()
seis_p=pressures
seis_rho, seis_vp, seis_vs, seis_vphi, seis_K, seis_G = \
    burnman.velocities_from_rock(rock, pressures, temperature, \
                                 burnman.averaging_schemes.VoigtReussHill())
depths= burnman.depths_for_rock(rock, pressures, temperature, \
                                 burnman.averaging_schemes.VoigtReussHill())

'''

print "preparations done"




# Priors on unknown parameters:
fraction_pv = pymc.Uniform('fraction_pv', 0.0, 1.0)
fe_pv = pymc.Uniform('fe_pv', 0.0, 0.5)
fe_pc= pymc.Uniform('fe_pc', 0.0, 0.5)
T0 = pymc.Normal('T0',mu=1500.,tau=1./500.)
def calc_all_velocities(fraction_pv, fe_pv, fe_pc,T0):
        method = 'slb3' #slb3|slb2|mgd3|mgd2
        pv=minerals.other.mg_fe_perovskite(fe_pv)
        pc=minerals.other.ferropericlase(fe_pc)
        rock = burnman.Composite( [fraction_pv, 1.0-fraction_pv],[ pv,pc ] )
        temperature = burnman.geotherm.adiabatic(seis_p,T0, rock)


        mat_rho, mat_vp, mat_vs, mat_vphi, mat_K, mat_G = burnman.velocities_from_rock(rock,seis_p, temperature)
        return mat_vp, mat_vs, mat_rho, mat_vphi, mat_K, mat_G

def nrmse(funca,funcb):
    """
    Normalized root mean square error for one profile
    :type funca: list of arrays of float
    :param funca: array calculated values
    :type funcb: list of arrays of float
    :param funcb: array of values (observed or calculated) to compare to

    :returns: RMS error
    :rtype: array of floats

    """
    diff=np.array(funca-funcb)
    diff=diff*diff
    rmse=np.sqrt(np.sum(diff)/len(diff))
    nrmse=rmse/(np.max(funcb)-np.min(funcb))
    return nrmse

def error(fraction_pv, fe_pv, fe_pc,T0):
    if True:
        if fraction_pv>0. and fraction_pv<1.0 and fe_pv>0. and fe_pv<1.0 and fe_pc>0. and fe_pc<1.0 and T0>300.:
            mat_vp, mat_vs, mat_rho, mat_vphi, mat_K, mat_G = calc_all_velocities(fraction_pv, fe_pv, fe_pc,T0)

            misfit = nrmse(mat_rho,seis_rho)+nrmse(mat_K,seis_K)+nrmse(mat_G,seis_G)
            return misfit
        else:
            return 1.e30

    #except:
    #    return 1e30

# rewrite as function of f_pv and fe_content
@pymc.deterministic
def calc_misfit(fraction_pv=fraction_pv, fe_pv=fe_pv, fe_pc=fe_pc,T0=T0):
    return error(fraction_pv,fe_pv,fe_pc,T0)


#sig = 1e-2
#misfit = pymc.Normal('d',mu=theta,tau=1.0/(sig*sig),value=0.,observed=True,trace=True)
#sig = pymc.Uniform('sig', 0.0, 100.0, value=1.)
sig=1.e-2 # Some sorts of error
obs = pymc.Normal('d',mu=calc_misfit,tau=1.0/(sig*sig),value=0,observed=True,trace=True)
model = [fraction_pv, fe_pv, fe_pc, T0,obs]
things = ['fraction_pv', 'fe_pv','fe_pc','T0']

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
    #pymc.MAP(model).fit() # Find minimum to start search from
    S = pymc.MCMC(model, db='pickle', dbname=dbname)
    S.sample(iter=100, burn=0, thin=1)
    S.db.close()
    whattodo="continue"

if whattodo=="continue":
    n_runs = 100
    for l in range(0,n_runs):
        db = pymc.database.pickle.load(dbname)
        print "*** run=%d/%d, # samples: %d" % (l, n_runs, db.trace('fraction_pv').stats()['n'] )
        S = pymc.MCMC(model, db=db)
        #S.sample(iter=100, burn=10, thin=1)
        S.sample(iter=10, burn=0, thin=1) # Search space for 100000 acceptable steps, forget first 1000 and save every 10.
        S.db.close()

if whattodo=="plot":
	files=sys.argv[2:]
	print "files:",files

        b=0#10000 # burn number
        i=1

        for t in things:
            print t
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
            print trace[0:100]
            for bin in [10,20,50,100]:
                hist,bin_edges=np.histogram(trace,bins=bin)
                a=np.argmax(hist)
                print "maxlike = ", bin_edges[a], bin_edges[a+1], (bin_edges[a]+bin_edges[a+1])/2.0


            (mu, sigma) = norm.fit(np.array(trace))
            print "mu, sigma: %e %e" % (mu, sigma)
            plt.subplot(2,(len(things)+1)/2,i)
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

    print "#samples: %d" % db.trace('fraction_pv').stats()['n']

    pymc.raftery_lewis(S, q=0.025, r=0.01)

    b = 1
    t = 1

    scores = pymc.geweke(S, intervals=20)
    print scores
    pymc.Matplot.trace(db.trace('deviance',chain=None).gettrace(burn=1000,thin=t,chain=None),'deviance',rows=2,columns=4,num=1)


    pymc.Matplot.trace(db.trace('fraction_pv',chain=None).gettrace(thin=t,chain=None),'fraction_pv',rows=2,columns=4,num=2)
    pymc.Matplot.histogram(np.array(db.trace('fraction_pv',chain=None).gettrace(burn=b,chain=None)),'fraction_pv',rows=2,columns=4,num=6)

    pymc.Matplot.trace(db.trace('fe_pv',chain=None).gettrace(thin=t,chain=None),'fe_pv',rows=2,columns=4,num=3)
    pymc.Matplot.histogram(np.array(db.trace('fe_pv',chain=None).gettrace(burn=b,chain=None)),'fe_pv',rows=2,columns=4,num=7)

    pymc.Matplot.trace(db.trace('fe_pc',chain=None).gettrace(thin=t,chain=None),'fe_pc',rows=2,columns=4,num=4)
    pymc.Matplot.histogram(np.array(db.trace('fe_pc',chain=None).gettrace(burn=b,chain=None)),'fe_pc',rows=2,columns=4,num=8)



    plt.show()


if whattodo=="show":
	values = [float(i) for i in sys.argv[2:]]
	#mat_vp, mat_vs, mat_rho = calc_velocities(values[0], values[1], values[2])
        mat_vp,mat_vs, mat_rho, mat_vphi=calc_all_velocities(values[0], values[1], values[2])
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
	[error(0.5,0.1,0.2) for i in range(0,1000)]

if whattodo=="profile":
	#just run normally
        print error(0.5,0.1,0.2)
	cProfile.run("[error(0.5,0.1,0.2) for i in range(0,1000)]")