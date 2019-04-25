import numpy as np
import scipy as sci
from scipy import integrate
import pylab
from scipy.integrate import odeint
from scipy import interpolate
from scipy.interpolate import CubicSpline
from scipy import misc

#Critical Temp=0.6175

T = 0.5955
k = 401
kb = 1.1
dofb = 3
nb = 1
kf = 1
doff = 4
nf = 1
cf = 13.94
cb = 16*cf
lam = 1/8
ms = np.sqrt(lam)
x = np.linspace(0, 2, k)
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
    return (((-((ms**2)/2)*x**2) + ((lam/4)*(x**4))) + (A*(dofb*((kb*x)**4))*(np.log(((kb*x)**2)/(Q**2))-1.5))\
           + (((dm2b / 2) * x ** 2) + ((dlb / 4) * x ** 4)))


def vfloop(x):
    return ((-((ms**2)/2)*x**2) + ((lam/4)*(x**4))) + (-A*((doff*((mf*x)**4))*(np.log(((mf*x)**2)/(Q**2))-1.5)))\
           + (((dm2f / 2) * x ** 2) + ((dlf / 4) * x ** 4))


def vtloop(x):
    return vbloop(x) + vfloop(x) - v0(x)


#Temperature corrections for potential

def intb(r, t):
    return (r**2) * np.log(1-np.exp(-np.sqrt((r**2) + (((kb*t)/T)**2))))


def intf(r, t):
    return (r**2) * np.log(1+np.exp(-np.sqrt((r**2) + (((kf*t)/T)**2))))


def solb(t):
    actb, errb = sci.integrate.quad(intb, 0.0, np.inf, args=(t))
    return actb*(dofb * ((T**4)/(2*np.pi**2)))


def solf(t):
    actf, errf = sci.integrate.quad(intf, 0.0, np.inf, args=(t))
    return actf*(-doff * ((T**4)/(2*np.pi**2)))


B = solf(x[1]) + solb(x[1])



def vtotal(t):
    return vtloop(t) + solb(t) + solf(t) - B


#Loop to plot vtotal:

vtotplot = []

for i in range(0, k):
    vtotplot.append(vtotal(x[i]))

#pylab.plot(x, vtotplot)
#pylab.show()


def dvtotal(t):
    return (vtotal(t+h) - vtotal(t-h)) / (2*h)

DV = []
VT = []
for i in range(0,k):
    DV.append(dvtotal(x[i]))
    VT.append(vtotal(x[i]))

#pylab.plot(x, DV, label='Diff')
#pylab.plot(x, VT, label='Potential')
#pylab.legend()
#pylab.show()

def du_dt(u, r):
    return [u[1], -2/(r+0.00001)*u[1] + dvtotal(u[0])]


u0 = [0.692, 0.0001]
xs = np.linspace(0, 100, k)
us = odeint(du_dt, u0, xs)
ys = us[:,0]
ysp = us[:,1]


#pylab.plot(xs, ys)
#pylab.show()

#Loop to plot vtotal(r):

vtotplotr = []

for i in range(0, k):
    vtotplotr.append(vtotal(ys[i]))

#Radial integral:


def radint(r):
    return r**2


#pylab.plot(xs, (vtotplotr + ysp**2/2)*radint(xs)*4*np.pi)
#pylab.show()


#Integrating the above graph to get s3/T:

s3int = []

for i in range(k):
    s3int.append((vtotplotr + 0.5 * ysp**2)*radint(xs))

s3 = []
s3 = np.trapz(s3int, xs)

#print(4*np.pi*s3/T)

#When T=0.61, s3/T=116.6 where u0=0.609772 and xs=75
#When T=0.6105, s3/T=130.08 where u0=0.60808 and xs=80.
#When T=0.61075, s3/T=138.02 where u0=0.6071396 and xs=82
#When T=0.611, s3/T=146.55 where u0=0.606105 and xs=85.
#When T=0.6125, s3/T=186.97 where u0=0.598555 and runs to xs=100.


p = [0.59, 0.595, 0.59525, 0.5955, 0.59575, 0.596, 0.6]
q = [63, 123.42, 128.15, 133.19, 138.45, 144.26, 337.3]

curve = CubicSpline(p, q)
xs1 = np.linspace(0.59, 0.6, k)
pylab.plot(xs1, curve(xs1))
pylab.show()
#pylab.plot(p, q)
#pylab.show()
nuctemp = sci.misc.derivative(curve, 0.5956, dx=0.001)
print(nuctemp)
print(curve(0.5956))

