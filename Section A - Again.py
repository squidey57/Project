import numpy as np
import scipy as sci
import pylab
from scipy import optimize
from scipy import integrate

k = 600
lam = 1/8
ms = np.sqrt(lam)
x = np.linspace(0, 2, k)
Q = 1
T = 0.6
A = 1/(64*np.pi**2)
mx = 1
dofb = 3
mf = 1
doff = 4
h = 0.01
xd = np.array([1, (1 + h), (1 - h)])

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


def v0(x):
    return (-((ms**2)/2)*x**2) + ((lam/4)*(x**4))


def vbloop(x):
    return(((-((ms**2)/2)*x**2) + ((lam/4)*(x**4))) + (A*(dofb*((mx*x)**4))*(np.log(((mx*x)**2)/(Q**2))-1.5)) \
           + (((dm2b / 2) * x ** 2) + ((dlb / 4) * x ** 4)))


def vfloop(x):
    return ((-((ms**2)/2)*x**2) + ((lam/4)*(x**4))) + (-A*((doff*((mf*x)**4))*(np.log(((mf*x)**2)/(Q**2))-1.5)))\
           + (((dm2f / 2) * x ** 2) + ((dlf / 4) * x ** 4))


def vtloop(x):
    return vbloop(x) + vfloop(x) - v0(x)


def intb(r):
    return (r**2) * np.log(1-np.exp(-np.sqrt((r**2) + ((mx*x[i])/T)**2)))

def intf(r):
    return (r**2) * np.log(1+np.exp(-np.sqrt((r**2) + ((mf*x[i])/T)**2)))

i3b = []
i3f = []
uv3b = np.zeros(k)
uv3f = np.zeros(k)
eb = np.zeros(k)
ef = np.zeros(k)
v3b = []
v3f = []
v3 = []
for i in range(0, k):
    i3b.append(sci.integrate.quad(intb,0,np.inf))
    i3f.append(sci.integrate.quad(intf,0,np.inf))
    uv3b[i], eb[i] = i3b[i]
    uv3f[i], ef[i] = i3f[i]
    v3b.append(dofb * ((T**4)/(2*np.pi**2)) * uv3b[i])
    v3f.append(-doff * ((T**4)/(2*np.pi**2)) * uv3f[i])
    v3.append(v3b[i] + v3f[i])

B = vtloop(x[1]) + v3[1]


def vtotal(x):
    return vtloop(x) + v3 - B

#pylab.plot(vtloop(x), label='vtloop')
#pylab.plot(vbloop(x), label='vbloop')
#pylab.plot(vfloop(x), label='vfloop')
pylab.plot(v0(x), label='v0')
pylab.plot(vtotal(x), label='Vtotal')
pylab.legend()
pylab.show()