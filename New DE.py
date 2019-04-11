import numpy as np
import scipy as sci
from scipy import integrate
import pylab
from scipy.integrate import odeint

T = 1
k = 401
kb = 1
dofb = 3
nb = 1
kf = 1
doff = 4
nf = 1
cf = 13.94
cb = 16*cf
lam = 1/8
ms = np.sqrt(lam)
x = np.linspace(-2, 2, k)
Q = 1
A = 1/(64*np.pi**2)
mx = 1
mf = 1
h = 0.01
xd = np.array([1, (1 + h), (1 - h)])


#Calculations for the quantum corrections
uv2b = []
uv2f = []

for i in range(0, 3):
    uv2b.append(A * dofb * ((mx * xd[i]) ** 4) * (np.log(((mx * xd[i]) ** 2) / Q ** 2) - 1.5))
    uv2f.append(A * -doff * ((mf * xd[i]) ** 4) * (np.log(((mf * xd[i]) ** 2) / Q ** 2) - 1.5))

v2x0b1 = uv2b[0]
v2xpb1 = uv2b[1]
v2xnb1 = uv2b[2]
dfv1b1 = ((v2xpb1 - v2xnb1) / (2 * h))
d2fv1b1 = ((v2xpb1 - (2 * v2x0b1) + v2xnb1) / (h ** 2))
dlb = 0.5 * (dfv1b1 - d2fv1b1)
dm2b = 0.5 * (1 / xd[0]) * ((d2fv1b1 * xd[0]) - (3 * dfv1b1))

v2x0f1 = uv2f[0]
v2xpf1 = uv2f[1]
v2xnf1 = uv2f[2]
dfv1f1 = (v2xpf1 - v2xnf1) / (2 * h)
d2fv1f1 = (v2xpf1 - (2 * v2x0f1) + v2xnf1) / (h ** 2)
dlf = 0.5 * (dfv1f1 - d2fv1f1)
dm2f = 0.5 * (1 / xd[0]) * ((d2fv1f1 * xd[0]) - (3 * dfv1f1))


#V-tree, pertubative and counter term corrections


def v0(x):
    return (-((ms**2)/2)*x**2) + ((lam/4)*(x**4))


def vbloop(x):
    return (((-((ms**2)/2)*x**2) + ((lam/4)*(x**4))) + (A*(dofb*((kb*x)**4))*(np.log(((kb*x)**2)/(Q**2))-1.5)) \
           + (((dm2b / 2) * x ** 2) + ((dlb / 4) * x ** 4)))


def vfloop(x):
    return ((-((ms**2)/2)*x**2) + ((lam/4)*(x**4))) + (-A*((doff*((mf*x)**4))*(np.log(((mf*x)**2)/(Q**2))-1.5)))\
           + (((dm2f / 2) * x ** 2) + ((dlf / 4) * x ** 4))


def vtloop(x):
    return vbloop(x) + vfloop(x) - v0(x)


#Temperature corrections for potential

def intb(r, t):
    return (r**2) * np.log(1-np.exp(-np.sqrt((r**2) + ((kb*t)/T)**2)))


def intf(r, t):
    return (r**2) * np.log(1+np.exp(-np.sqrt((r**2) + ((kf*t)/T)**2)))


def solb(t):
    actb, errb = sci.integrate.quad(intb, 0, np.inf, args=(t))
    return actb*(dofb * ((T**4)/(2*np.pi**2)))


def solf(t):
    actf, errf = sci.integrate.quad(intf, 0, np.inf, args=(t))
    return actf*(-doff * ((T**4)/(2*np.pi**2)))


B = solf(x[1]) + solb(x[1])


def vtotal(t):
    return vtloop(t) + solb(t) + solf(t) - B


def dvtotal(t):
    return (vtotal(t)[i+1] - vtotal(t)[i-1]) / (2*h)


def du_dr(u, r):
    return [u[1], (-2/r+0.0001)*u[1] + ]