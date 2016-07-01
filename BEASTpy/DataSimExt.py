### Data Generating Process
### Extended to include more complicated
### 1 juillet 2016
### Jeremy L Hour

import numpy as np
import math

def ClassicDataSim(n=2000,p=50,Ry=.5,Rd=.2,rho=.5):
  ### Covariate correlation coefficients
  Sigma = np.zeros(p*p)
  Sigma.shape = (p,p)
  
  for k in range(0,p):
    for j in range(0,p):
      Sigma[k,j] = rho**abs(k-j)

  ### Treatment effect
  a = 0
  
  ### Treatment variable coefficient
  gamma = np.zeros(p)
    
  for j in range(1,int(p/2)):
    gamma[j] = 1*(-1)**(j) / j**2
    
  ### Outcome equation coefficients
  b = gamma
    
  for j in range(int(p/2)+1,p):
    b[j] = (-1)**(j+1) / (p-j+1)**2
  
  ### Adjustment to match R.squared
  c = np.sqrt((1/gamma.dot(Sigma).dot(gamma))*(Rd/(1-Rd)))
  gamma = c*gamma
  
  c = math.sqrt((1/b.dot(Sigma).dot(b))*(Ry/(1-Ry)))
  b = c*b

  
  # Simulate covariates
  X = np.random.multivariate_normal(np.zeros(p), Sigma, n)
  
  # Simulate treatment
  d = np.random.uniform(size=n) < 1/(1+np.exp(-X.dot(gamma)))
  d = d.astype(int)
  
  # Simulate outcome
  y = a*d + X.dot(b) + np.random.normal(0,1,n)

  # Add the intercept
  X = np.c_[ np.ones((n,1)), X ]
  
  return X, y, d, b, gamma
  
  
def AwkwardDataSim(n=2000,p=50,Ry=.5,Rd=.2,rho=.5):
  ### Covariate correlation coefficients
  Sigma = np.zeros(p*p)
  Sigma.shape = (p,p)
  
  for k in range(0,p):
    for j in range(0,p):
      Sigma[k,j] = rho**abs(k-j)

  ### Treatment effect
  a = 0
  
  ### Treatment variable coefficient
  gamma = np.zeros(p)
    
  for j in range(1,int(p/2)):
    gamma[j] = 1*(-1)**(j) / j**2
    
  ### Outcome equation coefficients
  b = gamma
    
  for j in range(int(p/2)+1,p):
    b[j] = (-1)**(j+1) / (p-j+1)**2
  
  ### Adjustment to match R.squared
  c = np.sqrt((1/gamma.dot(Sigma).dot(gamma))*(Rd/(1-Rd)))
  gamma = c*gamma
  
  c = math.sqrt((1/b.dot(Sigma).dot(b))*(Ry/(1-Ry)))
  b = c*b

  
  # Simulate covariates
  X = np.random.multivariate_normal(np.zeros(p), Sigma, n)
  
  # Simulate treatment
  d = np.random.uniform(size=n) < 1/(1+np.exp(-X.dot(gamma)))
  d = d.astype(int)
  
  # Simulate outcome
  y = a*d + perturbation(X.dot(b)) + np.random.normal(0,1,n)

  # Add the intercept
  X = np.c_[ np.ones((n,1)), X ]
  
  return X, y, d, b, gamma 

def perturbation(x):
    return x * (1+.5*np.sin(x))