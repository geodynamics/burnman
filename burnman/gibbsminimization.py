# BurnMan - a lower mantle toolkit
# Copyright (C) 2012-2014, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

# NB, the dictionary of compositions -> array of compositions
# for endmembers and bulk composition should be done 
# during initialisation, not in the function itself

from scipy.optimize import fsolve

def null(A, eps=1e-15):
    u, s, vh = np.linalg.svd(A)
    s=np.append(s, [0. for i in range(len(vh)-len(s))])
    null_mask = (s <= eps)
    null_space = np.compress(null_mask, vh, axis=0)
    return np.transpose(null_space)

def equilibriumassemblage(bulk_composition, assemblage, pressure, temperature):
    # We can do this in one of two ways:
    # 1) Solve the equilibrium condition 
    # 2) Minimise the Gibbs Free Energy
    # These methods are mathematically, but not computationally equivalent

    # Method 1
    # Set up array of endmember compositions
    # Set up array of bulk composition
    # Find set of independent equations governing the equilibrium
    null=null(mbrcomp)
    mbrcomp=mbrcomp.T    
    bulk_composition=check_bulk(mbrcomp, comp)

# Create list of unknowns and their starting guesses
gvars=['P (kbar)', 'T (K)']
guesses=[250, 1000, 1.0] # First *mineral* proportion starts at 1.0

gvars.extend(['p(' + i + ')' for i in minlist])
guesses.extend(zeros(len(minlist)-1)) # All others start at 0.0

for i in range(len(lookup)):
    lookup[i]=[lookup[i][0]+len(guesses),lookup[i][1]+len(guesses)]

if len(gueslist) != 0:
    gvars.extend(zip(*gueslist)[0])
    guesses.extend(map(float,zip(*gueslist)[1]))

# Fix T, P, proportions or mixing parameters
fvars=[]
print '\nSelect two variables to fix'
for i in range(len(gvars)):
    print i, gvars[i]

for i in range(2):
    varinp=raw_input("Variable index and value (or min,max,nsteps): ")
    var=varinp.split()
    if len(var) == 2:
        fvars.append([int(var[0]), float(var[1]), float(var[1]), 1])
    elif len(var) == 4:
        fvars.append([int(var[0]), float(var[1]), float(var[2]), int(var[3])])
    else:
        print 'Incorrect number of lines'
        sys.exit()

print gvars[fvars[0][0]], fvars[0][1], '<= x <=', fvars[0][2], 'nsteps=', fvars[0][3]
print gvars[fvars[1][0]], fvars[1][1], '<= x <=', fvars[1][2], 'nsteps=', fvars[1][3]

# Set up problem and solve it with fsolve
print "\n", ' '.join(x.rjust(10) for x in gvars)
fixed_vars=[]
for i in linspace(fvars[0][1],fvars[0][2],fvars[0][3]):
    guesses[fvars[0][0]]=i
    if fvars[0][0] == 2:
        guesses[3]=1.0
    for j in linspace(fvars[1][1],fvars[1][2],fvars[1][3]):
        guesses[fvars[1][0]]=j
        fixed_vars=[[fvars[0][0], i],[fvars[1][0], j]]
        if fvars[1][0] == 2:
            guesses[3]=1.0
        elif fvars[0][0] == 2 and fvars[1][0] == 3:
            guesses[4]=1.0
        #guesses=[273.5461, 1385.0, 1.0, 0.0, 0.0, 0.0, 0.001306, 0.05165, 0.2963, 0.9044, 0.2589,0.9736,0.001794]
        #guesses=[253.7421,1219.776,1.0 ,  0.0 ,  0.0, 0.0 , 0.0 ,  0.2236  ,  0.9035 , 0.005607 , 0.001415,  0.007105  ,  0.2297, 0.1989 ]
        #print guesses
        soln=fsolve(set_eqns,guesses,args=(comp, minlist, mbrlist, mbrcomp, mbr, model, null, lookup, fixed_vars), full_output=1, xtol=1e-10)
        if soln[2]==1:
            print " ".join(str(("%10.5f" % x)) for x in soln[0])
            guesses=soln[0]
        else:
            print " ".join(str(("%10.5f" % x)) for x in soln[0]), soln[3]



def set_eqns(arg, comp, minlist, mbrlist, mbrcomp, mbr, model, null, lookup, fixed_vars):
    ''' 
    Equations to solve are:
    1) Complete set of independent equations
    2) Mass balance
    3) Zero mode 
    
    Unknowns are:
    1) Mineral proportions
    2) Solid solution compositions
    3) P and/or T
    
    arg[] is the list of variables for which to solve the equation
    comp is the elemental composition in mol %
    mbr is the complete set of minerals
    '''
    eqns=[]
    G=[]
    P=float(arg[0])
    T=float(arg[1]) # arg[0] and arg[1] are P, T
    mbrprp=zeros(len(mbrlist))
    for i in range(len(mbrlist)):
        if (mbrlist[i][1]==""): # Endmember, contained in "min"
            mbrprp[i]=arg[mbrlist[i][3]+2]
            G.append(gibbs(mbr[i],P,T))
        else: # Part of a solid solution
            n_mod=mbrlist[i][1]
            n_mbr=mbrlist[i][2]
            # print lookup[n_mod][0], lookup[n_mod][1]
            for j in range(lookup[n_mod][0],lookup[n_mod][1],1):
                if arg[j] > 1:
                    arg[j]=1.0
                if arg[j]<0:
                    arg[j]=0.0

            slnmbrprp=getprpns(model[n_mod], arg[lookup[n_mod][0]:lookup[n_mod][1]])[n_mbr]
            mbrprp[i]=arg[mbrlist[i][3]+2]*slnmbrprp
            if (mbr[i].H == "fictive"):
                Gpart=0.0
                for j in range(len(mbr[i].dpndnt)):
                    f=Fraction(model[mbrlist[i][1]].make[mbrlist[i][2]][2*j+2])
                    Gpart=Gpart + f*gibbs(mbr[i].dpndnt[j],P,T)
            else:
                Gpart=gibbs(mbr[i],P,T)
                

            # Add DQF
            if (model[mbrlist[i][1]].dqfs[mbrlist[i][2]] != []):
                dqf=model[mbrlist[i][1]].dqfs[mbrlist[i][2]]
                Gpart=Gpart + dqf[0] + dqf[1]*P + dqf[2]*T

            # Add RTlna
            if (slnmbrprp>1e-10):
                Gpart=Gpart + RTlna(model[n_mod], arg[lookup[n_mod][0]:lookup[n_mod][1]], n_mbr, P,T)

            G.append(Gpart)

    for i in range(len(null)): # Add independent equations
        eqns.append(dot(G,null[i]))

    # Fix two variables 
    eqns.append(arg[fixed_vars[0][0]]-fixed_vars[0][1])
    eqns.append(arg[fixed_vars[1][0]]-fixed_vars[1][1])
    
    # Add mineral composition constraints
    eqns.append(100.*(sum(arg[2:2+len(minlist)])-1.0)) # Add constraint that *mineral* proportions must sum to unity

    j=0
    for i in range(len(comp)):
        if (comp[i] != 0.):
            if j==0:
                factor=dot(mbrcomp[i],mbrprp)/comp[i]
            if j>0:
                if len(eqns)<len(arg):
                    #print elementorder[i], mbrcomp[i], mbrprp, comp[i], factor
                    eqns.append(10.*(dot(mbrcomp[i],mbrprp)-factor*comp[i]))
            j=j+1

    return eqns



def check_bulk(mbrcomp, comp):
    # Check that the bulk composition can be satisfied by a mixture of the prescribed minerals
    x,resid,rank,s = lstsq(mbrcomp,comp)
    resid = dot(mbrcomp,x) - comp
    tol=1e-10

    # Test cvxpy's ability to find a suitable range of endmembers
    #y = variable(len(x),1)
    #p = program(minimize(norm2(matrix(mbrcomp)*y-matrix(comp).T)),[equals(sum(y),1.0),geq(y,0.0)])
    #p.solve(quiet=True)
    #print y.value
    #print matrix(mbrcomp)*y.value-matrix(comp).T

    if (any(abs(resid[i]) > tol for i in range(len(resid)))):
        print "Oh no! Bulk composition outside the compositional range of the chosen minerals ( Maximum residual:", max(resid),  ")"

        tol=1e-3
        if (any(abs(resid[i]) > tol for i in range(len(resid)))):
            print "Residuals too high to recast composition. Program exiting"
            sys.exit()
        else:
            print "Residuals low enough to recast composition. Doing this now..."
            return dot(mbrcomp,x)

    return comp
