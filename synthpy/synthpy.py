# Synthetic Control functions definiion

import numpy as np
from cvxopt import matrix, solvers
from scipy.optimize import minimize

### Solution for a given V
def wsol(X0,X1,V):
    # For a given V, return the SC weights
    n = X0.shape[1]
    
    P = 2*(X0.T).dot(V.dot(X0))
    P = matrix(P.as_matrix())
    q = -2*X0.T.dot(V.dot(X1))
    q = matrix(q.as_matrix())
    
    ## Define constraints for cvx
    # Sum to one
    A = matrix(1.0,(1,n))
    b = matrix(1.0)
    # Between zero and one
    G = matrix( np.concatenate((np.identity(n),-np.identity(n))) )
    h = matrix( np.concatenate((np.ones(n),np.zeros(n))) )

    # Solving
    sol = solvers.qp(P,q,G,h,A,b)
    return sol['x']

### Function returning MSPE
def funcwsol(C,X0,X1,Z0,Z1):
    # Here C should be a vector of n-1
    Vw = np.insert(C,0,1)
    sol = wsol(X0,X1,np.diag(Vw))
    dis = Z1.as_matrix() - Z0.dot(sol)
    dis = dis.as_matrix()
    f = dis.T.dot(dis)
    return f
       
### Function to perform optimization over V
def synth(d,X,V=None,Z=None,OptimV=False):
    if V==None:
        V = np.identity(X.shape[1])
        
    X1 = X[d==1].T; X0 = X[d==0].T;
    sol = wsol(X0,X1,V)
    if OptimV:
        Z1 = Z[d==1].T; Z0 = Z[d==0].T
        p = X.shape[1]
        sol = minimize(funcwsol,np.ones(p-1),args=(X0,X1,Z0,Z1), method='nelder-mead', 
                       options={'xtol': 1e-8, 'disp': True})

    return sol
    
### L2 penalty for finding w solution    
def wsoll2(X0,X1,V,pen=0.0):
    # For a given V, return the SC weights
    n = X0.shape[1]
    
    # plus a L2 penalty term
    dis = X0 - np.kron(np.ones([1,n]), X1)
    
    P = 2*(X0.T).dot(V.dot(X0)) + 2*pen*(dis.T).dot(V.dot(dis))
    P = matrix(P.as_matrix())
    q = -2*X0.T.dot(V.dot(X1))
    q = matrix(q.as_matrix())
    
    ## Define constraints for cvx
    # Sum to one
    A = matrix(1.0,(1,n))
    b = matrix(1.0)
    # Between zero and one
    G = matrix( np.concatenate((np.identity(n),-np.identity(n))) )
    h = matrix( np.concatenate((np.ones(n),np.zeros(n))) )

    # Solving
    sol = solvers.qp(P,q,G,h,A,b)
    return sol['x']    
 
### Solution with l1 penalty
def wsol(X0,X1,V):
    # For a given V, return the SC weights
    n = X0.shape[1]
    
    P = 2*(X0.T).dot(V.dot(X0))
    P = matrix(P.as_matrix())
    q = -2*X0.T.dot(V.dot(X1))
    q = matrix(q.as_matrix())
    
    ## Define constraints for cvx
    # Sum to one
    A = matrix(1.0,(1,n))
    b = matrix(1.0)
    # Between zero and one
    G = matrix( np.concatenate((np.identity(n),-np.identity(n))) )
    h = matrix( np.concatenate((np.ones(n),np.zeros(n))) )

    # Solving
    sol = solvers.qp(P,q,G,h,A,b)
    return sol['x']   