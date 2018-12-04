import numpy as np
import pylab
import scipy as sci
from scipy import integrate
from scipy import misc
k = 1200
x = np.linspace(0, 1, k)
l = 1 / 8
m = np.sqrt(l)
kf = 2
dof = 4
kb = 1
dob = 3
A = 1 / (64 * (np.pi ** 2))
Q = 1


def v0(x):
    return -((m ** 2) / 2) * x ** 2 + (l / 4) * x ** 4


def v1f(x):
    return -A * dof * ((kf * x) ** 4) * ((np.log((kf * x) ** 2) / (Q ** 2)) - 1.5)


def v1b(x):
    return A * dob * ((kb * x) ** 4) * ((np.log((kb * x) ** 2) / (Q ** 2)) - 1.5)


##################
# V2 DIFFERENTIAL
##################

dfv1b1 = sci.misc.derivative(v1b,1,dx=10**-5, n =1)
d2fv1b1 = sci.misc.derivative(v1b,1,dx = 10**-5, n = 2)
dlb = 0.5 * (dfv1b1 - d2fv1b1)
dm2b = (0.5 * d2fv1b1 ) - (3 * dfv1b1)

dfv1f1 = sci.misc.derivative(v1f,1,dx=10**-5, n =1)
d2fv1f1 = sci.misc.derivative(v1f,1,dx = 10**-5, n = 2)
dlf = 0.5 * (dfv1f1 - d2fv1f1)
dm2f = (0.5 *d2fv1f1 ) - (3 * dfv1f1)
#############################################################
# V3 INTEGRATIONS#
T = 0.5248

i3b = []
i3f = []
uv3b = np.zeros(k)
uv3f = np.zeros(k)
eb = np.zeros(k)
ef = np.zeros(k)
v3b = []
v3f = []
v3 = []


def intb(r):
    return (r ** 2) * np.log(1 - np.exp(-np.sqrt((r ** 2) + ((kb * x[i]) / T) ** 2)))


def intf(r):
    return (r ** 2) * np.log(1 + np.exp(-np.sqrt((r ** 2) + ((kf * x[i]) / T) ** 2)))

for i in range(0, k):
    i3b.append(sci.integrate.quad(intb, 0, np.inf))
    i3f.append(sci.integrate.quad(intf, 0, np.inf))
    uv3b[i], eb[i] = i3b[i]
    uv3f[i], ef[i] = i3f[i]
    v3b.append((dob * ((T ** 4) / (2 * np.pi ** 2))) * uv3b[i])
    v3f.append((-dof * ((T ** 4) / (2 * np.pi ** 2))) * uv3f[i])
    v3.append(v3b[i] + v3f[i])


#################################################################
def v2f(x):
    return (dm2f / 2) * x ** 2 + (dlf / 4) * x ** 4


def v2b(x):
    return (dm2b / 2) * x ** 2 + (dlb / 4) * x ** 4


def vfloop(x):
    return v0(x) + v1f(x) + v2f(x)


def vbloop(x):
    return v0(x) + v1b(x) + v2b(x)


def vloop(x):
    return v0(x) + v1f(x) + v2f(x) + v1b(x) + v2b(x)

def v3f(t):
    return t*-dof * ((T ** 4) / (2 * np.pi ** 2)) * uv3f

def v3b(t):
    return t*dob * ((T ** 4) / (2 * np.pi ** 2)) * uv3f

def vtotal(x):
    return vloop(x) + v3 - B


def y(x):
    return x*0


B = vloop(x[1]) + v3[1]
#pylab.plot(y(x), label = 'V0')
#pylab.plot(vfloop(x), label = 'Fermion Loop')
#pylab.plot(vbloop(x), label = 'Boson Loop')
#pylab.plot(vloop(x), label = 'All Loop')
#pylab.plot(vtotal(x), label='Temperature')
#pylab.ylim(-0.00005, 0.00005)
#pylab.xlim(0, 30)
#pylab.legend()
#pylab.show()

vt1 = np.array([0.907, (684/1200)/0.656, (532/1200)/0.65, (296/1200)/0.6115, (123/1200)/0.5248])
k = np.array([1, 1.25,1.5, 1.75, 2])

pylab.plot(k,vt1)
pylab.show()


