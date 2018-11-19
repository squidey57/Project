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
    return -(A * dof[0] * (x ** 4) * (np.log((x ** 2 / Q ** 2) - 1.5)))
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
        v2f.append(A * -(dof[0] * ((x[n] * xd[i]) ** 4) * (np.log(((x * xd[i]) ** 2 / Q ** 2) - 1.5))))
        v2ba.append(A * dof[1] * ((mxa[n] * xd[i]) ** 4) * (np.log((((mxa[n] * xd[i]) ** 2) / (Q ** 2)) - 1.5)))
        v2bb.append(A * dof[2] * ((mxb[n] * xd[i]) ** 4) * (np.log((((mxb[n] * xd[i]) ** 2) / (Q ** 2)) - 1.5)))
        v2bc.append(A * dof[3] * ((mxc[n] * xd[i]) ** 4) * (np.log((((mxc[n] * xd[i]) ** 2) / (Q ** 2)) - 1.5)))
        v2bd.append(A * dof[4] * ((mxd[n] * xd[i]) ** 4) * (np.log((((mxd[n] * xd[i]) ** 2) / (Q ** 2)) - 1.5)))
        v2be.append(A * dof[5] * ((mxe[n] * xd[i]) ** 4) * (np.log((((mxe[n] * xd[i]) ** 2) / (Q ** 2)) - 1.5)))
#For boson A:
v2x0ba = v2ba[0]+v2f[0]
v2xpba = v2ba[1]+v2f[1]
v2xnba = v2ba[2]+v2f[2]

dfv2a = (v2xpba-v2x0ba)/h
d2fv2a = ((v2xpba - (2*v2x0ba) + v2xnba)/(h**2))

#dl and dm calculation:
dlba = 0.5*(((1/(xd[0]**3))*dfv2a) - ((1/(xd[0]**2))*d2fv2a))
dm2a = -d2fv2a - (3*dlba*xd[0]**2)

#For boson B:
v2x0bb = v2bb[0]+v2f[0]
v2xpbb = v2bb[1]+v2f[1]
v2xnbb = v2bb[2]+v2f[2]

dfv2b = (v2xpbb-v2x0bb)/h
d2fv2b = ((v2xpbb - (2*v2x0bb) + v2xnbb)/(h**2))

#dl and dm calculation:
dlbb = 0.5*(((1/(xd[0]**3))*dfv2b) - ((1/(xd[0]**2))*d2fv2b))
dm2b = -d2fv2b - (3*dlbb*xd[0]**2)

#For boson C:
v2x0bc = v2bc[0]+v2f[0]
v2xpbc = v2bc[1]+v2f[1]
v2xnbc = v2bc[2]+v2f[2]

dfv2c = (v2xpbc-v2x0bc)/h
d2fv2c = ((v2xpbc - (2*v2x0bc) + v2xnbc)/(h**2))

#dl and dm calculation:
dlbc = 0.5*(((1/(xd[0]**3))*dfv2c) - ((1/(xd[0]**2))*d2fv2c))
dm2c = -d2fv2c - (3*dlbc*xd[0]**2)

#For boson D:
v2x0bd = v2bd[0]+v2f[0]
v2xpbd = v2bd[1]+v2f[1]
v2xnbd = v2bd[2]+v2f[2]

dfv2d = (v2xpbd-v2x0bd)/h
d2fv2d = ((v2xpbd - (2*v2x0bd) + v2xnbd)/(h**2))

#dl and dm calculation:
dlbd = 0.5*(((1/(xd[0]**3))*dfv2d) - ((1/(xd[0]**2))*d2fv2d))
dm2d = -d2fv2d - (3*dlbd*xd[0]**2)

#For boson E:
v2x0be = v2be[0]+v2f[0]
v2xpbe = v2be[1]+v2f[1]
v2xnbe = v2be[2]+v2f[2]

dfv2e = (v2xpbe-v2x0be)/h
d2fv2e = ((v2xpbe - (2*v2x0be) + v2xnbe)/(h**2))

#dl and dm calculation:
dlbe = 0.5*(((1/(xd[0]**3))*dfv2e) - ((1/(xd[0]**2))*d2fv2e))
dm2e = -d2fv2e - (3*dlbe*xd[0]**2)

#Now to insert these new found dl and dm into the V2(counter term):

v2ba1 = []
v2bb1 = []
v2bc1 = []
v2bd1 = []
v2be1 = []

for i in range(0,k):
    v2ba1.append((((dm2a ** 2) / 2) * (x[i] ** 2)) + (((dlba ** 4) / 4) * (x[i] ** 4)))
    v2bb1.append((((dm2b ** 2) / 2) * (x[i] ** 2)) + (((dlbb ** 4) / 4) * (x[i] ** 4)))
    v2bc1.append((((dm2c ** 2) / 2) * (x[i] ** 2)) + (((dlbc ** 4) / 4) * (x[i] ** 4)))
    v2bd1.append((((dm2d ** 2) / 2) * (x[i] ** 2)) + (((dlbd ** 4) / 4) * (x[i] ** 4)))
    v2be1.append((((dm2e ** 2) / 2) * (x[i] ** 2)) + (((dlbe ** 4) / 4) * (x[i] ** 4)))


#Now to make the Vloop for the 5 bosons:
vloopa = []
vloopb = []
vloopc = []
vloopd = []
vloope = []

for i in range(0,k):
    vloopa.append(v0[i] + v1ba[i] + v2ba1[i])
    vloopb.append(v0[i] + v1bb[i] + v2bb1[i])
    vloopc.append(v0[i] + v1bc[i] + v2bc1[i])
    vloopd.append(v0[i] + v1bd[i] + v2bd1[i])
    vloope.append(v0[i] + v1be[i] + v2be1[i])


