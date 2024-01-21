import numpy as np

######Lagrange######
def lagrangeInterpolation(x, fx, target):
    Lx = np.zeros_like(target)
    for k in range(target.shape[0]):
        h = 0
        for i in range(x.shape[0]):
            za = 1
            for j in range(x.shape[0]):
                if i != j:
                    za = za*(target[k]-x[j])/(x[i]-x[j])
            za = za * fx[i]
            h += za
        Lx[k] = h
    return Lx

######Neville######
def nevilleInterpolation(x, fx, target):
    n = x.shape[0]
    res = np.zeros_like(target)
    for iter in range(target.shape[0]):
        Q = [[0] * n for _ in range(n)]
        for i in range(n):
            Q[i][0] = fx[i]
        for k in range(1, n):
            for i in range(n - k):
                Q[i][k] = ((target[iter] - x[i + k]) * Q[i][k - 1] + (x[i] - target[iter]) * Q[i + 1][k - 1]) / (x[i] - x[i + k])
        res[iter]=Q[0][n - 1]
    return res

######Hermite######
def dl(i, xi):
    result = 0.0
    for j in range(len(xi)):
        if j!=i:
            result += 1/(xi[i]-xi[j])
    return result

def l(i, xi, x):
    deno = 1.0
    nu = 1.0
    for j in range(len(xi)):
        if j!= i:
            deno *= (xi[i]-xi[j])
            nu *= (x-xi[j])
    return nu/deno

def hermiteInterpolation(x, fx, dyi, target):
    result = 0.0
    for i in range(len(x)):
        result += (fx[i]+(target-x[i])*(dyi[i]-2*fx[i]*dl(i, x)))*((l(i,x,target))**2)
    return result

######Newton######
def newtonInterpolation(x, fx, target):
    divdiff = getdivideddiff(x, fx)
    hatf = np.ones(len(target))*fx[0]
    for n in range(len(target)):
        t0 = 1
        for k in range(len(x)-1):
            # Newton
            t0 = (target[n]-x[k])*t0
            hatf[n] += divdiff[k]*t0
    return hatf

def getdivideddiff(x, fx):
    divdiff = np.zeros(len(x)-1)
    t = fx
    for k in range(len(x)-1):
        num0 = t        
        num1 = np.concatenate((np.ones(1),num0[0:-1]))    
        k0 = k+1
        den = np.concatenate((np.ones(k0),x[0:-k0]))
        t = (num0-num1)/(x-den)
        divdiff[k] = t[k0]
    return divdiff
