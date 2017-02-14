### Calibration vs. Synthetic Control
### Started while at Harvard
### Edited: 22 november 2016
### Jeremy L Hour

import os
import time
import numpy as np
import pandas as pd
from scipy.linalg import inv, solve, det
import scipy
import scipy.stats
import math
from cvxopt import matrix, solvers
from scipy.optimize import minimize
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

# Load functions
os.chdir("R:/Simulations/BEAST/synthpy")
from synthpy import wsol, funcwsol, synth
os.chdir("R:/Simulations/BEAST/BEASTpy")
from calibpy import gamma, gammagrad
from logitpy import logitloss, logitlossgrad
from DataSimExt import AwkwardDataSim
from ImmunizedATT import ImmunizedATT

# Monte Carlo simulation to evaluate calibration versus synthetic control

### Monte Carlo Simulation
R = 10000; p=5; n=100
V = np.identity(p)
val_b =[]; val_s = []; cf=[]
start_time = time.time()

for r in range(0,R):
    X, y, d, b0, g0 = AwkwardDataSim(n,p,Ry=.5,Rd=.2,rho=.5)
    #1. balancing
    res = minimize(gamma,np.zeros(p+1),args=(d,X), method='bfgs', jac=gammagrad,
              options={'gtol': 1e-8, 'disp': True})
    the = ImmunizedATT(y,d,X,res['x'], Immunity=False)
    val_b.append(the[0])

    #2. synthetic control
    M = pd.DataFrame(np.delete(X[d==1].mean(axis=0),0).T)
    X0 = pd.DataFrame(scipy.delete(X[d==0],0,1))
    sol = np.array(wsol(X0.T,M,V))
    the_s = np.mean(y[d==1]) - y[d==0].dot(sol)
    val_s.append(the_s[0])
    cf.append(y[d==0].dot(sol))

print("--- %s seconds ---" % (time.time() - start_time))
print(np.mean(val_b)); print(np.std(val_b))
print(np.mean(val_s)); print(np.std(val_s))

### Testing for normality
cf = np.array(cf)
test = scipy.stats.mstats.normaltest(cf)
print(test.pvalue)

test = scipy.stats.mstats.normaltest(val_s)
print(test.pvalue)

# check
W = np.exp(X[d==0].dot(res['x']))/sum(d) # compute calibration weights
BCc = np.squeeze(M) - X0.T.dot(W)
BCs = np.squeeze(M) - np.squeeze(X0.T.dot(sol))

# Plot of synthetic weights vs calibration
plt.scatter(W,sol)
plt.plot([0, .2], [0, .2], color = 'red', linestyle = 'solid')
plt.xlabel('Calibration weight')
plt.ylabel('Synthetic weight')
plt.title(r'$\mathrm{Calibration\ v\ Synthetic\ weights}$')
plt.show()

# Get corresponding beta for synthetic control
ly = np.log(sol*sum(d))
beta_SC = inv(X[d==0].T.dot(X[d==0])).dot(X[d==0].T.dot(ly))

### Distrbution of SC estimate
n, bins, patches = plt.hist(val_s, 60, normed=1, facecolor='red', alpha=0.5)
y = mlab.normpdf( bins, 0, np.std(val_s))
l = plt.plot(bins, y, 'r--', linewidth=2)

n, bins, patches = plt.hist(val_b, 60, normed=1, facecolor='blue', alpha=0.5)

#plot
plt.xlabel('ATT')
plt.ylabel('Probability')
plt.title(r'$\mathrm{Histogram\ of\ Synthetic\ Control\ (red)\ and\ Calibration\ (blue)\ ATT:}$')
plt.grid(True)

plt.show()

### Distrbution of SC counterfactual
n, bins, patches = plt.hist(cf, 60, normed=1, facecolor='red', alpha=0.5)
y = mlab.normpdf( bins, np.mean(cf), np.std(cf))
l = plt.plot(bins, y, 'r--', linewidth=2)
plt.xlabel('Y0')
plt.ylabel('Probability')
plt.title(r'$\mathrm{Histogram\ of\ Synthetic\ Control\ Outcome:}$')
plt.grid(True)
plt.show()

# regression plot
m, b = np.polyfit(val_b, val_s, 1)
plt.plot(val_b, val_s, '.')
plt.plot(val_b, m*val_b + b, '-')