### calibpy

import numpy as np

### Calibration objective and gradient
def gamma(beta,d,X):
    f = (1-d)*np.exp(X.dot(beta)) - d*(X.dot(beta))
    return np.mean(f)

def gammagrad(beta,d,X):
    g = (((1-d)*np.exp(X.dot(beta)) - d).T).dot(X)/X.shape[0]
    return g