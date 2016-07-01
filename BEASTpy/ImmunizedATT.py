### Immunized ATT

import numpy as np

def ImmunizedATT(y,d,X,beta,mu=None,Immunity=True):
    if mu==None:
        mu = np.zeros(X.shape[1])

    eps = y
    if Immunity:
        eps = y - X.dot(mu)

    pi = np.mean(d)
    theta = np.mean((d - (1-d)*np.exp(X.dot(beta)))*eps) / pi

    # Compute standard error
    psi = (d - (1-d)*np.exp(X.dot(beta)))*eps - d*theta
    VAR = np.mean(psi**2)/pi**2
    sigma = np.sqrt(VAR)/np.sqrt(X.shape[0])
      
    return theta, sigma