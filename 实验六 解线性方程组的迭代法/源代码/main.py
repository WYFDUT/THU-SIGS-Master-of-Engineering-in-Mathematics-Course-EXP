import time
import numpy as np 
import iterative as iter
import matplotlib.pyplot as plt
from scipy.sparse.linalg import bicgstab


if __name__ == "__main__":

        A, b = iter.generateMat(n=10000)
        #A, b = iter.generateSparseMat(n=10000)

        
        start = time.time()
        x, k = iter.jacobi(A,b,np.zeros_like(b),tol=1e-5)
        end = time.time()
        print(end-start,k)

        start = time.time()
        x,k = iter.gauss_seidel(A,b,np.zeros_like(b),tol=1e-5)
        end = time.time()
        print(end-start,k)

        # Compute Spectral Radius func
        #r = iter.cal_SpectralRadius(A, mode='gauss seidel')
        #print(r)

        start = time.time()
        x,k = iter.SOR(A,b,np.zeros_like(b), tol=1e-5, maxiter=10000, omega=1.2)
        end = time.time()
        print(end-start,k)

        start = time.time()
        x,k = iter.conjugate_gradient(A,b,np.zeros_like(b), tol=1e-5)
        end = time.time()
        print(end-start,k)

        start = time.time()
        solution, info = bicgstab(A, b, x0=np.zeros_like(b), tol=1e-5, maxiter=10000)
        end = time.time()
        print(solution, end-start)
        
                
        """
        e1, iteration1 = iter.jacobiWithError(A,b,np.zeros_like(b),tol=1e-5)
        e2, iteration2 = iter.gauss_seidelWithError(A,b,np.zeros_like(b),tol=1e-5)
        e3, iteration3 = iter.SORWithError(A,b,np.zeros_like(b), tol=1e-5, maxiter=10000, omega=1.3)
        plt.plot(iteration1, e1, c='blue')    
        plt.plot(iteration2, e2, c='orange')
        plt.plot(iteration3, e3, c='red')
        plt.legend(['Jacobi', 'Gauss Seidel', 'SOR'])

        plt.ylabel('inaccuracies', 
                        rotation=90,  
                )
        plt.xlabel('iterations', 
                )
        plt.title('Convergence curves of three iterative methods at n=100000')
        plt.show()
        """

        #start = time.time()
        #x, k = iter.jacobiSparse(A,b,np.zeros_like(b),tol=1e-5)
        #end = time.time()
        #print(end-start,k)

        #start = time.time()
        #x, k = iter.gauss_seidelSparse(A,b,np.zeros_like(b),tol=1e-5)
        #end = time.time()
        #print(end-start,k)

        #start = time.time()
        #x, k = iter.SORSparse(A,b,np.zeros_like(b), tol=1e-5, maxiter=10000, omega=1.3)
        #end = time.time()
        #print(end-start,k)





