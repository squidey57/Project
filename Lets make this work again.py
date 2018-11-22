import numpy as np
import scipy as sci
import pylab
from scipy import optimize

l = 1 / 8
m = np.sqrt(1 / 8)
x = np.linspace(0, 2, 400)
t = np.ones(400)
mx = 1
dox = 3
A = (1 / (64 * np.pi ** 2))
Q = 1
mf = 4
dof = 4

def v0(x):
    return (-((m ** 2) / 2) * x ** 2) + ((l / 4) * x ** 4)

def v1b(x):
    return A * dox * ((mx * x) ** 4) * (np.log(((mx * x ** 2) / Q ** 2)) - 1.5)

def v1f(x):
    return -A * dof * ((mf * x) ** 4) * (np.log(((mf * x ** 2) / Q ** 2)) - 1.5)

# pylab.plot(v0(x), label = 'VTree')
# pylab.legend()
# pylab.show()

# pylab.plot(v1b(x), label = 'V1B')
# pylab.legend()
# pylab.show()

# pylab.plot(v1f(x), label = 'V1F)
# pylab.legend()
# pylab.show()

h = 0.01
xd = np.array([1, (1 + h), (1 - h)])

uv2b = []
uv2f = []

for i in range(0, 3):
    uv2b.append(A * dox * ((mx * xd[i]) ** 4) * ((np.log(((mx * xd[i]) ** 2) / Q ** 2)) - 1.5))
    uv2f.append(-A * dof * ((mf * xd[i]) ** 4) * (np.log(((mf * xd[i] ** 2) / Q ** 2)) - 1.5))

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


def v2b(x):
    return ((dm2b / 2) * x ** 2) + ((dlb / 4) * x ** 4)


def v2f(x):
    return ((dm2f / 2) * x ** 2) + ((dlf / 4) * x ** 4)

#pylab.plot(v2b(x), label = 'V2B')
#pylab.plot(v2f(x), label = 'V2F')
#pylab.legend()
#pylab.show()

#pylab.plot(v0(x) + v1b(x)+v2b(x), label = 'Boson with Counter Term')
#pylab.plot(v0(x) + v1f(x)+v2f(x), label = 'Fermion with Counter Term')
#pylab.xlim(0,250)
#pylab.ylim(-0.5,0)
#pylab.legend()
#pylab.show()

T = 120
def intb(r):
    return (r**2) * np.log(1-np.exp(-np.sqrt((r**2) + ((mx*x[i])/T)**2)))

def intf(r):
    return (r**2) * np.log(1+np.exp(-np.sqrt((r**2) + ((mf*x[i])/T)**2)))

i3b = []
i3f = []
uv3b = np.zeros(400)
uv3f = np.zeros(400)
eb = np.zeros(400)
ef = np.zeros(400)

for i in range(0,400):
    i3b.append(sci.integrate.quad(intb,0,np.inf))
    i3f.append(sci.integrate.quad(intf,0,np.inf))
    uv3b[i], eb[i] = i3b[i]
    uv3f[i], ef[i] = i3f[i]




