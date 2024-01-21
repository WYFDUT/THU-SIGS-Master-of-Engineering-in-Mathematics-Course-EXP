import numpy as np


def trapz(f, a, b, N=50):
    x= np.linspace(a, b, N+1)
    y = f(x)
    y_right, y_left = y[1:], y[:-1]
    dx = (b-a) / N
    T = dx/2 * sum(y_right + y_left)
    return T


def simpson(f, a, b, N=50):
    dx, S = (b-a) / N, 0.0
    x = np.linspace(a, b, N+1)
    for i in range(N):
        S += (f(x[i]) + f(x[i+1]) + 4*f(x[i]+0.5*dx))
    return dx/6 * S


def Romberg(f, a, b, epsilon):
    h = b - a
    T = [[j for j in [0, 0, 0, 0]] for i in range(4)]
    i, j = 1, 1
    T[0][0] = h * (f(b) + f(a)) / 2
    T[1][0] = T[0][0] / 2 + f((b - a) / 2) * h / 2
    T[1][1] = (T[1][0] * 4 ** i - T[0][0]) / (4 ** i - 1)
    # # When row>4, until the accuracy reaches the given value
    while abs(T[i][i]-T[i-1][i-1]) > epsilon:
        i += 1
        if i == 4 : 
            break
        T[i][0] = trapz(f, a, b, N=2**(i))
        for j in range(1, i + 1):
            T[i][j] = (T[i][j - 1] * 4 ** i - T[i - 1][j - 1]) / (4 ** i - 1)
    if i < 4: 
        return T[i][j], i
    # When row>4, compute Romberg sequence 
    while abs(T[i-1][-1] - T[i-2][-1]) > epsilon:
        T.append([])
        T[i].append(trapz(f, a, b, N=2**(i)))
        for j in range(1, 5):
            T[i].append((T[i][j - 1] * 4 ** i - T[i - 1][j - 1]) / (4 ** i - 1))    
        i += 1
    return T[i-1][4], i-1


def simpleSimpson(f, a, b):
    return (b - a) / 6 * (f(a) + f(b) + 4 * f(a + 0.5 * (b - a)))

def Adaptive(f, a, b, epsilon):
    h = (a + b) / 2
    ST = simpleSimpson(f, a, b)
    SL = simpleSimpson(f, a, h)
    SR = simpleSimpson(f, h, b)
    if(abs(SL + SR - ST) < 15.0 * epsilon):
        return SL + SR + (SL + SR - ST) / 15.0
    return Adaptive(f, a, h, epsilon/2.0) + Adaptive(f, h, b, epsilon/2.0)
