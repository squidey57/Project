import numpy as np
import scipy as sci
import pylab
import matplotlib
from scipy import optimize
from scipy import integrate
matplotlib.use('Qt5Agg')
#############################################################################################
#############################################################################################
#DEFINING THE VARIABLES AND CONSTANTS
l = 1 / 8
m = np.sqrt(l)
k = 600
x = np.linspace(0, 1, k)
print(x)
t = np.ones(k)
mx = 1
dox = 3
A = (1 / (64 * np.pi ** 2))
Q = 1
mf = 0.5
dof = 4
h = 0.01
xd = np.array([1, (1 + h), (1 - h)])
#############################################################################################
#############################################################################################
#DIFFERENTIAL FOR V2
uv2b = []
uv2f = []

for i in range(0, 3):
    uv2b.append(A * dox * ((mx * xd[i]) ** 4) * (np.log(((mx * xd[i]) ** 2) / Q ** 2) - 1.5))
    uv2f.append(-A * dof * ((mf * xd[i]) ** 4) * (np.log(((mf * xd[i]) ** 2) / Q ** 2) - 1.5))

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
#############################################################################################

v0 = []
v1b = []
v1f = []
v2b = []
v2f = []
vtot = []
vnloop = []
i3b = []
i3f = []
uv3b = np.zeros(k)
uv3f = np.zeros(k)
eb = np.zeros(k)
ef = np.zeros(k)
v3b = []
v3f = []
v3 = []
vfloop = []
vbloop = []

#############################################################################################
#############################################################################################
#SETTING TEMPERATURE
T = 0.697
#############################################################################################
#############################################################################################

def intb(r):
    return (r ** 2) * np.log(1 - np.exp(-np.sqrt((r ** 2) + ((mx * x[i]) / T) ** 2)))


def intf(r):
    return (r ** 2) * np.log(1 + np.exp(-np.sqrt((r ** 2) + ((mf * x[i]) / T) ** 2)))



#############################################################################################
#############################################################################################
#Creating all the different potentials
# v0 = Tree Level Potential
# v1 = Perterbative Corrections to the Potential
# v2 = Counter Term Potential
# v3 = Thermal Correction
for i in range(0, k):
    i3b.append(sci.integrate.quad(intb, 0, np.inf))
    i3f.append(sci.integrate.quad(intf, 0, np.inf))
    uv3b[i], eb[i] = i3b[i]
    uv3f[i], ef[i] = i3f[i]
    v3b.append((dox * ((T ** 4) / (2 * np.pi ** 2))) * uv3b[i])
    v3f.append((-dof * ((T ** 4) / (2 * np.pi ** 2))) * uv3f[i])
    v3.append(v3b[i] + v3f[i])
    v0.append(-(((m ** 2) / 2) * x[i] ** 2) + ((l / 4) * x[i] ** 4))
    v1b.append((A * dox * ((mx * x[i]) ** 4) * ((np.log((mx * x[i] ** 2)/(Q**2))) - 1.5)))
    v1f.append((-A * dof * ((mf * x[i]) ** 4) * ((np.log((mf * x[i] ** 2)/(Q**2))) - 1.5)))
    v2b.append((((dm2b / 2) * x[i] ** 2) + ((dlb / 4) * x[i] ** 4)))
    v2f.append((((dm2f / 2) * x[i] ** 2) + ((dlf / 4) * x[i] ** 4)))
    vtot.append(v0[i] + v1f[i] + v1b[i] + v2f[i] + v2b[i] + v3[i])
    vnloop.append(v0[i] + v1f[i] + v1b[i] + v2f[i] + v2b[i])
    vbloop.append(v0[i] + v1b[i] + v2b[i])
    vfloop.append(v0[i] + v1f[i] + v2f[i])

def y(x):
    return 0*x
#############################################################################################
#############################################################################################
#PLOTTING GRAPHS


pylab.plot(y(x), label='V0')
#pylab.plot(vbloop, label = 'Vbloop')
#pylab.plot(vfloop, label = 'Vfloop')
pylab.plot(vtot - vtot[1], label='VTotal')
#pylab.ylim(-50, 250)
#pylab.xlim(0, 10)
pylab.legend()
pylab.show()

k = np.array([0.25, 0.5, 0.75, 1, 1.25, 1.5])
vt = np.array([(374.5/600)/0.7425, (363/600)/0.697])

pylab.plot(k[0:2], vt)
pylab.show()