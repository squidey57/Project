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
v0 = []
for i in range(0,k):
    v0.append(-(((ms**2)/2)*x[i]**2)+((l/4)*x[i]**4))

#V1 Calculations:

v1f = []
v1ba = []
v1bb = []
v1bc = []
v1bd = []
v1be = []

for i in range(0,k):
    v1f.append(-(A * dof[0] * (x[i] ** 4) * ((np.log((x[i] ** 2 / Q ** 2))) - 1.5)))
    v1ba.append(A * dof[1] * (x[i] ** 4) * (np.log(((x[i]) ** 2) / (Q ** 2)) - 1.5))
    v1bb.append(A * dof[1] * ((2*x[i])**4) * ((np.log(((2*x[i])**2 )/ (Q ** 2))) - 1.5))
    v1bc.append(A * dof[1] * ((3*x[i])**4) * ((np.log(((3*x[i])**2 )/ (Q ** 2))) - 1.5))
    v1bd.append(A * dof[1] * ((4*x[i])**4) * ((np.log(((4*x[i])**2 )/ (Q ** 2))) - 1.5))
    v1be.append(A * dof[1] * ((5*x[i])**4) * ((np.log(((5*x[i])**2 )/ (Q ** 2))) - 1.5))


#V2 calculations:
h = 0.001
xd = np.array([256,256+h,256-h])


v2f = []
v2ba = []
v2bb = []
v2bc = []
v2bd = []
v2be = []
for i in range(0,3):
    v2f.append(A * -(dof[0] * (( xd[i]) ** 4) * (np.log((( xd[i]) ** 2 / Q ** 2) - 1.5))))
    v2ba.append(A * dof[1] * (( xd[i]) ** 4) * (np.log(((( xd[i]) ** 2) / (Q ** 2)) - 1.5)))
    v2bb.append(A * dof[2] * ((2 * xd[i]) ** 4) * (np.log((((2* xd[i]) ** 2) / (Q ** 2)) - 1.5)))
    v2bc.append(A * dof[3] * ((3 * xd[i]) ** 4) * (np.log((((3* xd[i]) ** 2) / (Q ** 2)) - 1.5)))
    v2bd.append(A * dof[4] * ((4 * xd[i]) ** 4) * (np.log((((4*xd[i]) ** 2) / (Q ** 2)) - 1.5)))
    v2be.append(A * dof[5] * ((5 * xd[i]) ** 4) * (np.log((((5* xd[i]) ** 2) / (Q ** 2)) - 1.5)))
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
    vloopa.append(v0[i] + v1ba[i] + v2ba1[i] + v1f[i])
    vloopb.append(v0[i] + v1bb[i] + v2bb1[i] + v1f[i])
    vloopc.append(v0[i] + v1bc[i] + v2bc1[i] + v1f[i])
    vloopd.append(v0[i] + v1bd[i] + v2bd1[i] + v1f[i])
    vloope.append(v0[i] + v1be[i] + v2be1[i] + v1f[i])

##########################################################################
T = 249.809
#Ta = 249.809 , Va=124

def ifer(r):
    return (r**2)*np.log(1+np.exp(-np.sqrt((r**2)+((((x[i]))**2)/T**2))))
def iba(r):
    return (r**2)*np.log(1-np.exp(-np.sqrt((r**2)+((((x[i])))**2)/T**2)))
def ibb(r):
    return (r**2)*np.log(1-np.exp(-np.sqrt((r**2)+(((2*(x[i])))**2)/T**2)))
def ibc(r):
    return (r**2)*np.log(1-np.exp(-np.sqrt((r**2)+(((3*(x[i])))**2)/T**2)))
def ibd(r):
    return (r**2)*np.log(1-np.exp(-np.sqrt((r**2)+(((4*(x[i])))**2)/T**2)))
def ibe(r):
    return (r**2)*np.log(1-np.exp(-np.sqrt((r**2)+(((5*(x[i])))**2)/T**2)))


intifer = []
intiba = []
intibb = []
intibc = []
intibd = []
intibe = []


for i in range(0,k):
    intifer.append(sci.integrate.quad(ifer,0,np.inf))
    intiba.append(sci.integrate.quad(iba,0,np.inf))
    intibb.append(sci.integrate.quad(ibb, 0, np.inf))
    intibc.append(sci.integrate.quad(ibc, 0, np.inf))
    intibd.append(sci.integrate.quad(ibd, 0, np.inf))
    intibe.append(sci.integrate.quad(ibe, 0, np.inf))

uv3f = np.zeros(k)
uv3a = np.zeros(k)
uv3b = np.zeros(k)
uv3c = np.zeros(k)
uv3d = np.zeros(k)
uv3e = np.zeros(k)

ef = np.zeros(k)
ea = np.zeros(k)
eb = np.zeros(k)
ec = np.zeros(k)
ed = np.zeros(k)
ee = np.zeros(k)


for i in range(0,k):
    uv3f[i], ef[i] = intifer[i]
    uv3a[i], ea[i] = intiba[i]
    uv3b[i], eb[i] = intibb[i]
    uv3c[i], ec[i] = intibc[i]
    uv3d[i], ed[i] = intibd[i]
    uv3e[i], ee[i] = intibe[i]

v3f = dof[0]*((T**4)/(2*np.pi**2))*uv3f
v3a = dof[1]*((T**4)/(2*np.pi**2))*uv3a
v3b = dof[2]*((T**4)/(2*np.pi**2))*uv3b
v3c = dof[3]*((T**4)/(2*np.pi**2))*uv3c
v3d = dof[3]*((T**4)/(2*np.pi**2))*uv3d
v3e = dof[3]*((T**4)/(2*np.pi**2))*uv3e


v3a1 = []
v3b1 = []
v3c1 = []
v3d1 = []
v3e1 = []

for i in range(0,k):
    v3a1.append(v3f[i] + v3a[i])
    v3b1.append(v3f[i] + v3b[i])
    v3c1.append(v3f[i] + v3c[i])
    v3d1.append(v3f[i] + v3d[i])
    v3e1.append(v3f[i] + v3e[i])

vtota = []
vtotb = []
vtotc = []
vtotd = []
vtote = []


for i in range(0,k):
    vtota.append(v3a1[i] + vloopa[i])
    vtotb.append(v3b1[i] + vloopb[i])
    vtotc.append(v3c1[i] + vloopc[i])
    vtotd.append(v3d1[i] + vloopd[i])
    vtote.append(v3e1[i] + vloope[i])
print(min(vtota[30:k]-vtota[1]))
pylab.plot(vtota - vtota[1])
pylab.ylim(1,-1)
pylab.xlim(120,128)
pylab.show()


