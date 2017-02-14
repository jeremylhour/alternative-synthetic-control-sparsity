### Two ways of using SC with micro data

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
from synthpy import wsol, funcwsol, synth, wsoll2
os.chdir("R:/Simulations/BEAST/BEASTpy")
from calibpy import gamma, gammagrad
from logitpy import logitloss, logitlossgrad
from DataSimExt import AwkwardDataSim
from ImmunizedATT import ImmunizedATT

### Comparison of two ways
p=5; n=100
V = np.identity(p)
X, y, d, b0, g0 = AwkwardDataSim(n,p,Ry=.5,Rd=.2,rho=.5)
val_b =[]; val_s = []

#1. balancing
res = minimize(gamma,np.zeros(p+1),args=(d,X), method='bfgs', jac=gammagrad,
               options={'gtol': 1e-8, 'disp': True})
the = ImmunizedATT(y,d,X,res['x'], Immunity=False)

#2. synthetic control
X0 = pd.DataFrame(scipy.delete(X[d==0],0,1))
for i in range(0,n):
    if d[i]==1:
        M = pd.DataFrame(np.delete(X[i],0).T)
        sol = np.array(wsol(X0.T,M,V))
        the_s = np.mean(y[d==1]) - y[d==0].dot(sol)
        val_s.append(the_s[0])
    else:
        print("not treated")
        

plt.scatter(y[d==1], val_s)
plt.show()
# check
W = np.exp(X[d==0].dot(res['x']))/sum(d) # compute calibration weights
BCc = np.squeeze(M) - X0.T.dot(W)
BCs = np.squeeze(M) - np.squeeze(X0.T.dot(sol))



#### Testing the penalized version
p=5; n=100
V = np.identity(p)
X, y, d, b0, g0 = AwkwardDataSim(n,p,Ry=.5,Rd=.2,rho=.5)

M = pd.DataFrame(np.delete(X[d==1].mean(axis=0),0).T)
X0 = pd.DataFrame(scipy.delete(X[d==0],0,1))

R = 1000
penalty = np.linspace(0,10e7,R)
att = []

for val in penalty:
    sol = np.array(wsoll2(X0.T,M,V,pen=val))
    att.append(np.mean(y[d==1]) - y[d==0].dot(sol))

att = np.array(att)
### PLot solution as penalty varies
plt.axes().set_aspect('equal')
plt.xlabel('Penalty')
plt.ylabel('ATT')
plt.title('ATT as penalty varies')
 
# plot sample data
plt.plot(penalty,att,'b--',label='ATT')

