import numpy as np
import scipy as sci
from scipy import integrate
x =np.linspace(0,10,11)
b = np.array([0,1,2,3,4])
def v0(x):
    return x**2 + b[i]
int = []
for i in range(0,5):
    int.append(sci.integrate.quad(v0,0,5))

def v1(x):
    return x**2 + 4
print(int)
print(sci.integrate.quad(v1, 0, 5))

