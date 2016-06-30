### Python for small dimension case

### Define functions
def gamma(beta,d,X):
    f = (1-d)*np.exp(X.dot(beta)) - d*(X.dot(beta))
    return np.mean(f)

def gammagrad(beta,d,X):
    g = (((1-d)*np.exp(X.dot(beta)) - d).T).dot(X)/X.shape[0]
    return g

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

### Logit log-likelihood
def logitloss(beta,d,X):
    d_tilde = 2*d - 1
    f = np.log(1+np.exp(-d_tilde*(X.dot(beta))))
    return np.mean(f)

def logitlossgrad(beta,d,X):
    d_tilde = 2*d - 1
    f = -d_tilde / (1+np.exp(d_tilde*(X.dot(beta))))
    f = (f.T).dot(X)/X.shape[0]
    return f


### Start using them
X, y, d, b0, gamma0 = DataSim(n=3000,p=5,Ry=.5,Rd=.2,rho=.5)

from scipy.optimize import minimize
import time

# balancing
res = minimize(gamma,np.zeros(6),args=(d,X), method='nelder-mead',
              options={'xtol': 1e-8, 'disp': True})
b0 = np.insert(b0,0,0)

gamma(res['x'], d, X)
gamma(b0, d, X)

ImmunizedATT(y,d,X,res['x'], Immunity=False)

# logit
res = minimize(logitloss,np.zeros(6),args=(d,X), method='nelder-mead',
              options={'xtol': 1e-8, 'disp': True})
logitloss(res['x'], d, X)
logitloss(b0, d, X)

ImmunizedATT(y,d,X,res['x'], Immunity=False)

### Monte Carlo Simulation
R = 1000
val_b =[]; val_l = []
start_time = time.time()

for r in range(0,R):
    X, y, d, b0, g0 = DataSim(n=10000,p=10,Ry=.5,Rd=.2,rho=.5)
    #1. balancing
    res = minimize(gamma,np.zeros(11),args=(d,X), method='bfgs', jac=gammagrad,
              options={'gtol': 1e-8, 'disp': True})
    the = ImmunizedATT(y,d,X,res['x'], Immunity=False)
    val_b.append(the)

    #2. logit
    res = minimize(logitloss,np.zeros(11),args=(d,X), method='bfgs', jac=logitlossgrad,
              options={'gtol': 1e-8, 'disp': True})
    the = ImmunizedATT(y,d,X,res['x'], Immunity=False)
    val_l.append(the)

print("--- %s seconds ---" % (time.time() - start_time))
print(np.mean(val_b)); print(np.std(val_b))
print(np.mean(val_l)); print(np.std(val_l))


## Lasso example with cvxpy
from cvxpy import *
import matplotlib.pyplot as plt

# Problem data.
n = 15
m = 10
numpy.random.seed(1)
A = numpy.random.randn(n, m)
b = numpy.random.randn(n, 1)
# gamma must be positive due to DCP rules.
gamma = Parameter(sign="positive")

# Construct the problem.
x = Variable(m)
error = sum_squares(A*x - b)
obj = Minimize(error + gamma*norm(x, 1))
prob = Problem(obj)

# Construct a trade-off curve of ||Ax-b||^2 vs. ||x||_1
sq_penalty = []
l1_penalty = []
x_values = []
gamma_vals = numpy.logspace(-4, 6)
for val in gamma_vals:
    gamma.value = val
    prob.solve()
    # Use expr.value to get the numerical value of
    # an expression in the problem.
    sq_penalty.append(error.value)
    l1_penalty.append(norm(x, 1).value)
    x_values.append(x.value)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.figure(figsize=(6,10))

# Plot trade-off curve.
plt.subplot(211)
plt.plot(l1_penalty, sq_penalty)
plt.xlabel(r'\|x\|_1', fontsize=16)
plt.ylabel(r'\|Ax-b\|^2', fontsize=16)
plt.title('Trade-Off Curve for LASSO', fontsize=16)

# Plot entries of x vs. gamma.
plt.subplot(212)
for i in range(m):
    plt.plot(gamma_vals, [xi[i,0] for xi in x_values])
plt.xlabel(r'\gamma', fontsize=16)
plt.ylabel(r'x_{i}', fontsize=16)
plt.xscale('log')
plt.title(r'\text{Entries of x vs. }\gamma', fontsize=16)

plt.tight_layout()
plt.show()
