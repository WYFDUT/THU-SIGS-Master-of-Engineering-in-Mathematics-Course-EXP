import numpy as np
from scipy.sparse import diags


####### Mat Generate ########
def generateMat(n=100):
    array_a = np.int16(np.diag([3] * n))
    array = np.diag([-1] * (n-1))
    a, b = np.zeros((n-1)), np.zeros(n)
    array_b = np.insert(array, 0, values=a, axis=0)
    array_b = np.insert(array_b, (n-1), values=b, axis=1)
    array_c = np.insert(array, (n-1), values=a, axis=0)
    array_c = np.insert(array_c, 0, values=b, axis=1)
    array_d = np.fliplr(np.eye(n))*1/2
    matrix_A = array_a + array_b + array_c + array_d
    vector_b = np.ones(n) * 1.5
    if n%2 == 0:
        matrix_A[int(n/2-1), int(n/2)], matrix_A[int(n/2), int(n/2-1)] = -1, -1
        vector_b[0], vector_b[-1], vector_b[int(n/2-1)], vector_b[int(n/2)] = 2.5, 2.5, 1.0, 1.0
    else:
        matrix_A[int(n/2), int(n/2)] = 3
        vector_b[0], vector_b[-1], vector_b[int(n/2)] = 2.5, 2.5, 1.0
    return matrix_A, vector_b

def generateSparseMat(n):
    # Create diagonal matrix
    diag_a = np.ones(n) * 3
    diag_b = np.ones(n-1) * (-1)
    diagonals = [diag_a, diag_b, diag_b]
    # Diagonal matrix offsets
    offsets = [0, 1, -1]  # Adjusted offsets for sub-diagonal
    # Generate sparse matrix
    matrix_A = diags(diagonals, offsets, shape=(n, n), format='lil')
    for i in range(n):
        if matrix_A[i, n-1-i] != -1 and matrix_A[i, n-1-i] != 3:
            insertElement(matrix_A, i, n-1-i, 0.5)
    # Set values at specific positions
    # Create right-hand vector
    vector_b = np.ones(n) * 1.5
    if n % 2 == 0:
        vector_b[0], vector_b[-1], vector_b[int(n/2-1)], vector_b[int(n/2)] = 2.5, 2.5, 1.0, 1.0
    else:
        vector_b[0], vector_b[-1], vector_b[int(n/2)] = 2.5, 2.5, 1.0
    return matrix_A, vector_b

def insertElement(matrix, row, col, value):
    # Insert an element into the sparse matrix
    matrix[row, col] = value


####### Normal method ########
def jacobi(A, b, x0, tol=1e-6, maxiter=10000):
    n = len(A)
    x = x0.copy()
    for k in range(maxiter):
        x_old = x.copy()
        for i in range(n):
            x[i] = (b[i] - np.dot(A[i,:], x_old) + A[i,i] * x_old[i]) / A[i,i]
        if np.linalg.norm(x - x_old) < tol:
            break
    return x, k+1


def gauss_seidel(A, b, x0, tol=1e-3, maxiter=10000):
    n = len(A)
    x = x0.copy()
    for k in range(maxiter):
        x_old = x.copy()
        for i in range(n):
            x[i] = (b[i] - np.dot(A[i, :i], x[:i]) - np.dot(A[i, i+1:], x_old[i+1:])) / A[i, i]
        if np.linalg.norm(x - x_old) < tol:
            break
    return x, k+1


def SOR(A, b, x0, tol=1e-6, maxiter=10000, omega=1.2):
    n = len(A)
    x = x0.copy()
    for k in range(maxiter):
        for i in range(n):
            x[i] = (1 - omega) * x[i] + (omega / A[i, i]) * (
                        b[i] - np.dot(A[i, :i], x[:i]) - np.dot(A[i, i + 1:], x0[i + 1:]))
        if np.linalg.norm(x - x0) < tol:
            break
        x0 = x.copy()
    return x, k+1


def cal_SpectralRadius(A, mode='jacobi', omega=1.0):
    U = np.triu(A, k=1)
    L = np.tril(A, k=-1)
    D = A - L - U
    if mode == 'jacobi':
        B = np.dot(np.linalg.inv(D), (-(L+U)))
        spectral_radius = np.max(np.abs(np.linalg.eigvals(B)))
        return spectral_radius
    elif mode == 'gauss seidel':
        B = np.dot(np.linalg.inv(D+L), (-(U)))
        spectral_radius = np.max(np.abs(np.linalg.eigvals(B)))
        return spectral_radius
    elif mode == 'SOR':
        B = np.dot(np.linalg.inv(D+omega*L), ((1-omega)*D-omega*U))
        spectral_radius = np.max(np.abs(np.linalg.eigvals(B)))
        return spectral_radius


def jacobiWithError(A, b, x0, tol=1e-6, maxiter=10000):
    n = len(A)
    x = x0.copy()
    error, iter = [], []
    for k in range(maxiter):
        x_old = x.copy()
        for i in range(n):
            x[i] = (b[i] - np.dot(A[i,:], x_old) + A[i,i] * x_old[i]) / A[i,i]
        if np.linalg.norm(x - x_old) < tol:
            error.append(np.linalg.norm(x - x_old))
            iter.append(k+1)
            break
        else:
            error.append(np.linalg.norm(x - x_old))
            iter.append(k+1)
    return error, iter


def gauss_seidelWithError(A, b, x0, tol=1e-3, maxiter=10000):
    n = len(A)
    x = x0.copy()
    error, iter = [], []
    for k in range(maxiter):
        x_old = x.copy()
        for i in range(n):
            x[i] = (b[i] - np.dot(A[i, :i], x[:i]) - np.dot(A[i, i+1:], x_old[i+1:])) / A[i, i]
        if np.linalg.norm(x - x_old) < tol:
            error.append(np.linalg.norm(x - x_old))
            iter.append(k+1)
            break
        else:
            error.append(np.linalg.norm(x - x_old))
            iter.append(k+1)
    return error, iter


def SORWithError(A, b, x0, tol=1e-6, maxiter=10000, omega=1.2):
    n = len(A)
    x = x0.copy()
    error, iter = [], []
    for k in range(maxiter):
        for i in range(n):
            x[i] = (1 - omega) * x[i] + (omega / A[i, i]) * (
                        b[i] - np.dot(A[i, :i], x[:i]) - np.dot(A[i, i + 1:], x0[i + 1:]))
        if np.linalg.norm(x - x0) < tol:
            error.append(np.linalg.norm(x - x0))
            iter.append(k+1)
            break
        else:
            error.append(np.linalg.norm(x - x0))
            iter.append(k+1)
        x0 = x.copy()
    return error, iter


####### Sparse method ########
def jacobiSparse(A, b, x0, tol=1e-6, maxiter=10000):
    n = len(b)
    x = x0  
    for k in range(maxiter):
        x_new = np.zeros_like(x)
        for i in range(n):
            x_new[i] = (b[i] - A[i, :].dot(x) + A[i, i] * x[i]) / A[i, i]
        if np.linalg.norm(x_new - x) < tol:
            break
        x = x_new
    return x_new, k + 1

def gauss_seidelSparse(A, b, x0, tol=1e-3, maxiter=10000):
    n = len(b)
    x = x0  
    for k in range(maxiter):
        for i in range(n):
            x[i] = (b[i] - A[i, :i].dot(x[:i]) - A[i, i+1:].dot(x[i+1:])) / A[i, i]
        if np.linalg.norm(A.dot(x) - b) < tol:
            break
    return x, k + 1

def SORSparse(A, b, x0, tol=1e-6, maxiter=10000, omega=1.2):
    n = len(b)
    x = x0.copy()
    for k in range(maxiter):
        for i in range(n):
            x[i] = (1 - omega) * x[i] + (omega / A[i, i]) * (
                        b[i] - A[i, :i].dot(x[:i]) - A[i, i+1:].dot(x[i+1:]))
        if np.linalg.norm(x - x0) < tol:
            break
        x0 = x.copy()
    return x, k+1


###### Conjugate Gradient ######
def conjugate_gradient(A, b, x0, tol=1e-6, maxiter=10000):
    n = len(A)
    x = x0.copy()
    r = b - np.dot(A, x)
    p = r.copy()
    for k in range(maxiter):
        alpha = np.dot(r, r) / np.dot(p, np.dot(A, p))
        x = x + alpha * p
        r_new = r - alpha * np.dot(A, p)
        beta = np.dot(r_new, r_new) / np.dot(r, r)
        p = r_new + beta * p
        if np.linalg.norm(r_new) < tol:
            break
        r = r_new
    return x, k+1