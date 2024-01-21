import os
import numpy as np
import interpolation as Inter
import matplotlib.pyplot as plt


if __name__ == "__main__":
    x = np.array([1,2,3,4,5,6,7])
    fx = lambda x: 1/x
    #fx = np.array([0.368,0.135,0.050,0.018,0.007,0.002,0.001])
    fx_div = lambda x: -1/x**2
    deriv = fx_div(x)
    target = np.arange(0.05, 7.05, 0.05)
    y = Inter.newtonInterpolation(x, fx(x), target)

    plt.figure()
    plt.plot(target, fx(target))
    plt.plot(target, y, c='orange')
    for i in range(len(x)):    
        plt.plot(x[i], (fx(x[i])),'ro',markerfacecolor='none')
        #plt.plot(x[i], (fx[i]),'ro',markerfacecolor='none')
    plt.legend(['original','Lagrange','control point'])
    #plt.legend(['Newton','Lagrange','control point'])
    plt.title('Lagrange Interpolation')
    plt.show()