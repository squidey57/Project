import numpy as np
import scipy as sci
from scipy import integrate
from scipy import optimize
import pylab

A = (1/(64*np.pi**2))
k = 400
v = 256
x = np.linspace(0,k,k)
kf = 1
nf = kf
l = 1/8
ms = (np.sqrt(2*l)*v)
nx = np.array([1,2,3,4,5])
kx = np.array([1,2,3,4,5])
dof = np.array([4,3,6,9,12,15])
Q = 173



#V0 Calculations:
def v0(x):
    return -(((ms**2)/2)*x**2)+((l/4)*x**4)

#V1 Calculations:
mx = []
for n in range(0,5):
    for i in range(0,k):
        mx.append(kx[n] * x[i])

mxa = mx[0:400]
mxb = mx[400:800]
mxc = mx[800:1200]
mxd = mx[1200:1600]
mxe = mx[1600:2000]

def v1f(x):
    return A * dof[0] * (x ** 4) * (np.log((x ** 2 / Q ** 2) - 1.5))
def v1ba(x):
    return A * dof[1] * ((mxa) ** 4) * (np.log((((mxa) ** 2) / (Q ** 2)) - 1.5))
def v1bb(x):
    return A * dof[2] * ((mxb)**4) * (np.log((((mxb)**2 )/ (Q ** 2)) - 1.5))
def v1bc(x):
    return A * dof[3] * ((mxc)**4) * (np.log((((mxc)**2 )/ (Q ** 2)) - 1.5))
def v1bd(x):
    return A * dof[4] * ((mxd)**4) * (np.log((((mxd)**2 )/ (Q ** 2)) - 1.5))
def v1be(x):
    return A * dof[5] * ((mxe)**4) * (np.log((((mxe)**2 )/ (Q ** 2)) - 1.5))


#V2 calculations:
h = 0.001
xd = np.array([256,256+h,256-h])


v2f = []
v2ba = []
v2bb = []
v2bc = []
v2bd = []
v2be = []
for n in range(0,k):
    for i in range(0,3):
        v2f.append(A * dof[0] * ((x[n] * xd[i]) ** 4) * (np.log(((x * xd[i]) ** 2 / Q ** 2) - 1.5)))
        v2ba.append(A * dof[1] * ((mxa[n] * xd[i]) ** 4) * (np.log((((mxa[n] * xd[i]) ** 2) / (Q ** 2)) - 1.5)))
        v2bb.append(A * dof[2] * ((mxb[n] * xd[i]) ** 4) * (np.log((((mxb[n] * xd[i]) ** 2) / (Q ** 2)) - 1.5)))
        v2bc.append(A * dof[3] * ((mxc[n] * xd[i]) ** 4) * (np.log((((mxc[n] * xd[i]) ** 2) / (Q ** 2)) - 1.5)))
        v2bd.append(A * dof[4] * ((mxd[n] * xd[i]) ** 4) * (np.log((((mxd[n] * xd[i]) ** 2) / (Q ** 2)) - 1.5)))
        v2be.append(A * dof[5] * ((mxe[n] * xd[i]) ** 4) * (np.log((((mxe[n] * xd[i]) ** 2) / (Q ** 2)) - 1.5)))


