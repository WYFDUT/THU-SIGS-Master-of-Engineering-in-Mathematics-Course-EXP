import os 
import OLS
import numpy as np
import matplotlib.pyplot as plt


if __name__ =="__main__":
    t0=np.array([ 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55])
    y0=np.array([ 1.27, 2.16, 2.86, 3.44, 3.87, 4.15, 4.37, 4.51, 4.58, 4.02, 4.64])
    t=np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55])
    y=np.array([0, 1.27, 2.16, 2.86, 3.44, 3.87, 4.15, 4.37, 4.51, 4.58, 4.02, 4.64])
    target = np.arange(-5, 61, 1)
    res = OLS.expLeastSquare(t0, y0, target, 1)
    res2 = OLS.expLeastSquare(t0, y0, t, 1)
    print(np.mean(np.abs(res2-y)))

    plt.figure()
    plt.plot(target, res, c='orange')
    
    for i in range(len(t)):
        plt.arrow(t[i],res2[i],0,y[i]-res2[i],linewidth=1,color='g' ,head_length=0.05,head_width=0.5)
        plt.plot(t[i], res2[i],'yo')    
        plt.plot(t[i], y[i],'ro',markerfacecolor='none')
    plt.legend(['fitting curve', 'error', 'fitting point', 'sample point'])
    #plt.legend(['fitting point', 'sample point'])

    plt.xlabel('t(min)', 
                rotation=0
            )
    plt.ylabel('(x10-4)', 
                rotation=90,  
                horizontalalignment='right'
            )
    plt.title('Carbon content versus time')
    plt.show()
