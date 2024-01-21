import numpy as np


def ordinaryLeastSquare(x, y, target, power):
    A = np.ones(x.shape)
    res = np.zeros_like(target).astype("float64")
    for i in range(1, power+1):
        A=np.c_[A, x**i]
    A2 = A.T.dot(A)
    print(A2.shape)
    y2 = A.T.dot(y)
    param = np.linalg.inv(A2).dot(y2)
    for i in range(power+1):
        print(i)
        res += param[i]*target**i 
    print(param)
    return res

def expLeastSquare(x, y, target, power):
    y = np.log(y)
    x = np.log(x)
    A = np.ones(x.shape)
    res = np.zeros_like(target).astype("float64")
    A=np.c_[A, x]
    A2 = A.T.dot(A)
    y2 = A.T.dot(y)
    param = np.linalg.inv(A2).dot(y2)
    res += np.exp(param[0])*target**(param[1])
    return res
