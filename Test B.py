import numpy as np
import scipy as sci
import pylab
from scipy import integrate
from scipy import optimize

l = 1 / 8
j = 400


m = np.sqrt(2 * l)
A = (1 / (64 * np.pi ** 2))
Q = 173
x = np.linspace(0, 2, j)
z = 1
d = 0
nx = 1
kx = 1
mx = 1
dox = 3

nf = np.array([ 1, 2, 3, 4, 5])
kf = np.array([ 1, 2, 3, 4, 5])
mf = np.array([ 1, 2, 3, 4, 5])
dof = 4 * nf

v0 = []
for i in range(0, j):
    v0.append(-(((m ** 2) / 2) * x[i] ** 2) + ((l / 4) * x[i] ** 4))

v1x = []
v1fd = []


for i in range(0, j):
    v1x.append(A * dox * ((mx * x[i]) ** 4) * (np.log(((mx * x[i]) ** 2) / Q ** 2)) - 1.5)
    v1fd.append(-(A * dof[z] * ((mf[d] * x[i]) ** 4) * (np.log(((mf[d] * x[i]) ** 2) / Q ** 2)) - 1.5))


# Requires the calculation of a first and second derivative
h = 0.001
xd = np.array([246, 246 + h, 246 - h])

# It is possible to do the derivatives of V2 algebraically, but V1 is not as easy to calculate
# The masses will be slightly different for +h and -h for each case, so they must calculated seperately


# Now to calculate a +h and -h for the V1
uv2x = []
uv2f1 = []


for i in range(0, 3):
    uv2x.append(A * dox * ((mx * xd[i]) ** 4) * (np.log(((mx * xd[i]) ** 2) / Q ** 2)) - 1.5)
    uv2f1.append(-(A * dof[z] * ((mf[d] * xd[i]) ** 4) * (np.log(((mf[d] * xd[i]) ** 2) / Q ** 2)) - 1.5))


v2x0f1 = uv2x[0] + uv2f1[0]
v2xpf1 = uv2x[1] + uv2f1[1]
v2xnf1 = uv2x[2] + uv2f1[2]
dfv1f1 = (v2xpf1 - v2x0f1) / h
d2fv1f1 = (v2xpf1 - (2 * v2x0f1) + v2xnf1) / (h ** 2)
dlf1 = 0.5 * (((1 / (xd[0] ** 3)) * dfv1f1) - ((1 / (xd[0] ** 2)) * d2fv1f1))
dm2f1 = -d2fv1f1 - (3 * dlf1 * xd[0] ** 2)

v2fd = []

for i in range(0, j):
    v2fd.append(((dm2f1 / 2) * x[i] ** 2) + ((dlf1 / 4) * x[i] ** 4))


vloopfd = []


for i in range(0, j):
    vloopfd.append(v0[i] + v1x[i] + v1fd[i] + v2fd[i])

Tx = 1


def ix(r):
    return (r ** 2) * np.log(1 - np.exp(-np.sqrt((r ** 2) + (x[i] / Tx) ** 2)))


def if1(r):
    return (r ** 2) * np.log(1 + np.exp(-np.sqrt((r ** 2) + ((mf[0] * x[i]) / Tx) ** 2)))



intx = []
intfd = []

for i in range(0, j):
    intx.append(sci.integrate.quad(ix, 0, np.inf))
    intfd.append(sci.integrate.quad(if1, 0, np.inf))

print(intx)
uv3x = np.zeros(j)
uv3fd = np.zeros(j)

ex = np.zeros(j)
ef1 = np.zeros(j)

for i in range(0, j):
    uv3x[i], ex[i] = intx[i]
    uv3fd[i], ef1[i] = intfd[i]


v3x = dox * ((Tx ** 4) / (2 * np.pi ** 2)) * uv3x
nv3fd = -dof[z] * ((Tx ** 4) / (2 * np.pi ** 2)) * uv3fd

v3f1 = []
v3f2 = []
v3f3 = []
v3f4 = []
v3f5 = []

for i in range(0, j):
    v3f1.append(v3x[i] + nv3fd[i])

vtotfd = []


for i in range(0, j):
    vtotfd.append(vloopfd[i] + v3f1[i])

print(max(vtotfd[30:j]-vtotfd[1]))
#pylab.plot(vtotfd - vtotfd[1])
#pylab.ylim(-10,10)
#pylab.xlim(300,320)
#pylab.show()

pylab.plot(v0)
pylab.show()
#c1 = np.array([(199/212.79), (340/233.3873872), (317/330), (314/491.85), (313/712.3)])#
#c2 = np.array([(182/184.5), (320/224.23), (313/351.15235), (313/551.6569388), (314/818.4707081)])
#pylab.plot(kf[0:5],c1,  label = 'Nf = 1')
#pylab.legend()
#pylab.savefig('Nf = 1')
#pylab.show()