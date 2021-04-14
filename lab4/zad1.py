import numpy as np
from math import log

def f(x):
    return (1/(1+x**2))

def rectangle():
    length = 1
    h=0.1
    xs = np.linspace(0, 1, num=10, endpoint=False)
    Int = sum(map(lambda x: f(x)*h, xs[1::]))

    print("wartość: ", Int, "bład bezwględny: ", Int - np.pi/4)

def trapez():

    length = 1
    h = 0.1
    xs = np.linspace(0, 1, num=11, endpoint=True)
    Int = sum(map(lambda x: (f(x)+f(x-h))*h/2, xs[1::]))

    print("wartość: ", Int, "bład bezwględny: ", Int - np.pi/4)

def simpson():
    length = 1
    h = 0.1
    xs = np.linspace(0, 1, num=6, endpoint=True)
    Int = sum(map(lambda x: (f(x)+4*f(x-h)+f(x-2*h))*(h/3), xs[1::]))

    print("wartość: ", Int, "bład bezwględny: ", Int - np.pi/4)


rectangle()
print("---------")
trapez()
print("---------")
simpson()