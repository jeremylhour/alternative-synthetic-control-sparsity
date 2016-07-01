### logitpy

import numpy as np

def logitloss(beta,d,X):
    d_tilde = 2*d - 1
    f = np.log(1+np.exp(-d_tilde*(X.dot(beta))))
    return np.mean(f)

def logitlossgrad(beta,d,X):
    d_tilde = 2*d - 1
    f = -d_tilde / (1+np.exp(d_tilde*(X.dot(beta))))
    f = (f.T).dot(X)/X.shape[0]
    return f