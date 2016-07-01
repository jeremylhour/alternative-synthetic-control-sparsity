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
from synthpy import wsol, funcwsol, synth


### California Tobacco example
# Load synth dataset
data = pd.read_table("/Volumes/USB_KEY/BEAST/synthpy/calitobacco.txt")
data.head()
d = data['Treated']
X = data[["Income","RetailPrice", "Young", "BeerCons"
                  , "SmokingCons1970", "SmokingCons1971", "SmokingCons1972", "SmokingCons1973", "SmokingCons1974", "SmokingCons1975"
                  , "SmokingCons1980", "SmokingCons1988"]]
Z = data[["SmokingCons1970", "SmokingCons1971", "SmokingCons1972", "SmokingCons1973", "SmokingCons1974", "SmokingCons1975",
          "SmokingCons1976", "SmokingCons1977", "SmokingCons1978", "SmokingCons1979","SmokingCons1980"]]

V=np.diag([.1,.1,.1,.1,6,6,6,6,6,6,6,6])

sol = wsol(X[d==0].T,X[d==1].T,V)
print(X[d==1].T)
print(X[d==0].T.dot(sol))

sol = synth(d,X,V)

### Simulated data