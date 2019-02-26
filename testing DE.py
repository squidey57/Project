import numpy as np
import scipy as sci
import pylab
from scipy import optimize
from scipy import integrate
from scipy.integrate import odeint
from sympy import *

h = 0.01
div = []
x = np.linspace(0, 10, 100)


def f(x):
    return x**2
def a(x):
    return((f(x + h) - f(x - h)) / (2 * h))


pylab.plot(a(x))
pylab.plot(f(x))
pylab.show()