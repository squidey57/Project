import numpy as np
import scipy as sci
from scipy import integrate
k = 400
j = k+1
x = np.linspace(0,2,j)
kb = 1
kf = 1
dob = 3
dof = 4
T = 1

i3b = []
i3f = []
uv3b = np.zeros(j)
uv3f = np.zeros(j)
eb = np.zeros(j)
ef = np.zeros(j)


def mb(x):
    return kb*x

def mf(x):
    return kf*x

#print(mb(x))
def intb(r):
    return (r ** 2) * np.log(1 - np.exp(-np.sqrt((r ** 2) + (mb(x)[i] / T) ** 2)))

def intf(r):
    return (r ** 2) * np.log(1 + np.exp(-np.sqrt((r ** 2) + ((kf * x[i]) / T) ** 2)))

for i in range(0, j):
    i3b.append(sci.integrate.quad(intb, 0, np.inf))
    i3f.append(sci.integrate.quad(intf, 0, np.inf))
    uv3b[i], eb[i] = i3b[i]
    uv3f[i], ef[i] = i3f[i]

def v3f(t):
    return t*dob * ((T ** 4) / (2 * np.pi ** 2)) * uv3f

print(v3f(1))
