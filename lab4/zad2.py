import numpy as np
from scipy.integrate import quadrature

def f(x):
    return 1/(1+x**2)


def genLegendre(n):
    p = []
    p.append(np.poly1d([1]))
    p.append(np.poly1d([1,0]))
    x = np.poly1d([1,0])
    for i in range(1,n):
        p.append((2*i+1)/(i+1)*p[i]*x-(i/(i+1))*p[i-1])
    return p

def approxF(n):
    L = genLegendre(n)
    c = []
    fA = []
    fA.append(np.poly1d([0]))
    for i in range(0,n):
        c.append(quadrature((lambda x: f(x)*L[i](x)), -1,1)[0]/quadrature((lambda x: L[i](x)**2), -1, 1)[0])
        fA.append(c[i]*L[i])
    return sum(fA)

def Int(n):
    i = approxF(n).integ()
    return i(1) - i(-1)

print("Wartość obliczona: ", Int(8), "błąd bezwględny: ", Int(8) - np.pi/2)            
