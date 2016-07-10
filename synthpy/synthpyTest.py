### Synthpy test
### Jeremy L Hour
### 30 juin 2016

# Load modules
import os
import time
import numpy as np
import pandas as pd
from scipy.linalg import inv, solve, det
import scipy.stats
import math
from cvxopt import matrix, solvers
from scipy.optimize import minimize
import matplotlib.pyplot as plt

# Load functions
os.chdir("/Volumes/USB_KEY/BEAST/synthpy")
from synthpy import wsol, funcwsol, synth, wsoll2


### California Tobacco example
# Load synth dataset
data = pd.read_table("/Volumes/USB_KEY/BEAST/synthpy/calitobacco.txt")
data.head()
d = data['Treated']
X = data[["Income","RetailPrice", "Young", "BeerCons"
                  , "SmokingCons1970", "SmokingCons1971", "SmokingCons1972", "SmokingCons1973", "SmokingCons1974", "SmokingCons1975"
                  , "SmokingCons1980", "SmokingCons1988"]]

names = []
for i in range(1970,1981):
    names.append('SmokingCons' +str(i))
         
Z = data[names]

for i in range(1981,2001):
    names.append('SmokingCons' +str(i))
       
y = data[names]

### Data loading over
V=np.diag([.1,.1,.1,.1,6,6,6,6,6,6,6,6])

sol = wsol(X[d==0].T,X[d==1].T,V)
print(X[d==1].T)
print(X[d==0].T.dot(sol))

sol = synth(d,X,V)

### With L2 penalty term
R = 1000
penalty = np.linspace(0,1000000,R)
M = []

for val in penalty:
    sol = wsoll2(X[d==0].T,X[d==1].T,V,pen=val)
    a = y[d==1].as_matrix().T - y[d==0].T.dot(sol)
    M.append(np.array(a))
M = np.squeeze(M)

u_bound = M.max(axis=0)
l_bound = M.min(axis=0)

# set-up the plot
plt.axes().set_aspect('equal')
plt.xlabel('Year')
plt.ylabel('ATT')
plt.title('ATT as penalty varies')
 
# plot sample data
plt.plot(np.linspace(1970,2000,31),u_bound,'b--',label='Upper bound')
plt.plot(np.linspace(1970,2000,31),l_bound,'b--',label='Lower bound')
 
plt.show()