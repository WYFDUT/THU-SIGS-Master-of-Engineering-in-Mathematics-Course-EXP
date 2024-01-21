import time 
import numpy as np
from scipy.integrate import quad
import integration as Inte
import matplotlib.pyplot as plt


if __name__ == "__main__":
    time_ada, time_romberg = [], []

    f = lambda x: np.sqrt(4-(np.sin(x))**2)
    #f = lambda x: np.sin(x)/x
    #f = lambda x: np.exp(x)*np.sin(x)
    #f = lambda x: np.log(1+x)/(1+x**2)

    start = time.time()
    #a1 = Inte.trapz(f, 0, np.pi/4, N=80)
    #a1 = Inte.trapz(f, 1e-10, 1, N=80)
    #a1 = Inte.trapz(f, 1, 3, N=80)
    #a1 = Inte.trapz(f, 0, 1, N=80)
    end = time.time()
    print(end-start)

    start = time.time()
    #a2 = Inte.simpson(f, 0, np.pi/4, N=80)
    #a2 = Inte.simpson(f, 1e-10, 1, N=80)
    #a2 = Inte.simpson(f, 1, 3, N=80)
    #a2 = Inte.simpson(f, 0, 1, N=80)
    end = time.time()
    print(end-start)

    start = time.time()
    #a3 = Inte.Romberg(f, 0, np.pi/4, 1e-10)
    #a3 = Inte.Romberg(f, 1e-10, 1, 1e-10)
    #a3 = Inte.Romberg(f, 1, 3, 1e-10)
    #a3 = Inte.Romberg(f, 0, 1, 1e-10)
    end = time.time()
    print(end-start)

    for i in range(1,11):
        start = time.time()
        #a4 = Inte.Adaptive(f, 0, np.pi/4, 1e-5)
        #a4 = Inte.Adaptive(f, 1e-10, 1, 1e-5)
        a4 = Inte.Adaptive(f, 1, 3, 10**(-i))
        #a4 = Inte.Adaptive(f, 0, 1, 1e-5)
        end = time.time()
        time_ada.append(end-start)

        start = time.time()
        a3 = Inte.Romberg(f, 1, 3, 10**(-i))
        end = time.time()
        time_romberg.append(end-start)
        #print(end-start)

    #real = quad(f, 0, np.pi/4)
    #real = quad(f, 1e-10, 1)
    #real = quad(f, 1, 3)
    #real = quad(f, 0, 1)
    #print(a1, a2, a3, real, abs(a1-real[0]), abs(a2-real[0]), abs(a3[0]-real[0]), abs(a4-real[0]))
        
    # Plot Adaptive versus Romberg time cost contrast
    t0=np.array(['1e-1', '1e-2', '1e-3', '1e-4', '1e-5', '1e-6', '1e-7', '1e-8', '1e-9', '1e-10'])
    plt.figure()
    plt.plot(t0, time_ada, c='blue')    
    plt.plot(t0, time_romberg, c='orange')
    plt.legend(['Adaptive Simpson', 'Romberg'])
    #plt.legend(['fitting point', 'sample point'])

    plt.ylabel('(s)', 
                rotation=90,  
                horizontalalignment='right'
            )
    plt.title('Comparison of the running time of the two algorithms')
    plt.show()