import numpy as np
import scipy as sci
import pylab
from scipy import integrate
from scipy import optimize


l = 1/8
j = 600
v = 246
m = np.sqrt(2*l)*v
A = (1/(64*np.pi**2))
Q = 173
x = np.linspace(0,j,j)

nx = 1
kx = 1
mx = 1
dox = 3

nf = np.array([0, 1, 2, 3, 4, 5])
kf = np.array([0, 1, 2, 3, 4, 5])
mf = np.array([0, 1, 2, 3, 4, 5])
dof = 4 * nf


v0 = []
for i in range(0,j):
    v0.append(-(((m**2)/2)*x[i]**2) + ((l/4)*x[i]**4))

v1x = []
v1f1 = []
v1f2 = []
v1f3 = []
v1f4 = []
v1f5 = []

for m in range(1,6):
    for i in range(0,j):
        v1x.append(A * dox * ((mx*x[i])**4) * (np.log(((mx*x[i])**2)/Q**2))  - 1.5 )
        v1f1.append(-(A * dof[1] * ((mf[m] * x[i]) ** 4) * (np.log(((mf[m] * x[i]) ** 2) / Q ** 2)) - 1.5))
        v1f2.append(-(A * dof[2] * ((mf[m] * x[i]) ** 4) * (np.log(((mf[m] * x[i]) ** 2) / Q ** 2)) - 1.5))
        v1f3.append(-(A * dof[3] * ((mf[m] * x[i]) ** 4) * (np.log(((mf[m] * x[i]) ** 2) / Q ** 2)) - 1.5))
        v1f4.append(-(A * dof[4] * ((mf[m] * x[i]) ** 4) * (np.log(((mf[m] * x[i]) ** 2) / Q ** 2)) - 1.5))
        v1f5.append(-(A * dof[5] * ((mf[m] * x[i]) ** 4) * (np.log(((mf[m] * x[i]) ** 2) / Q ** 2)) - 1.5))


#Requires the calculation of a first and second derivative
h = 0.001
xd = np.array([246,246+h,246-h])

#It is possible to do the derivatives of V2 algebraically, but V1 is not as easy to calculate
#The masses will be slightly different for +h and -h for each case, so they must calculated seperately


#Now to calculate a +h and -h for the V1
uv2x = []
uv2f1 = []
uv2f2 = []
uv2f3 = []
uv2f4 = []
uv2f5 = []
for m in range(1,6):
    for i in range(0,3):
        uv2x.append(A * dox * ((mx*xd[i])**4) * (np.log(((mx*xd[i])**2)/Q**2))  - 1.5 )
        uv2f1.append(-(A * dof[1] * ((mf[m] * xd[i]) ** 4) * (np.log(((mf[m] * xd[i]) ** 2) / Q ** 2)) - 1.5))
        uv2f2.append(-(A * dof[2] * ((mf[m] * xd[i]) ** 4) * (np.log(((mf[m] * xd[i]) ** 2) / Q ** 2)) - 1.5))
        uv2f3.append(-(A * dof[3] * ((mf[m] * xd[i]) ** 4) * (np.log(((mf[m] * xd[i]) ** 2) / Q ** 2)) - 1.5))
        uv2f4.append(-(A * dof[4] * ((mf[m] * xd[i]) ** 4) * (np.log(((mf[m] * xd[i]) ** 2) / Q ** 2)) - 1.5))
        uv2f5.append(-(A * dof[5] * ((mf[m] * xd[i]) ** 4) * (np.log(((mf[m] * xd[i]) ** 2) / Q ** 2)) - 1.5))

v2x0f1 = uv2x[0]+uv2f1[0]
v2xpf1 = uv2x[1]+uv2f1[1]
v2xnf1 = uv2x[2]+uv2f1[2]
dfv1f1 = (v2xpf1-v2x0f1)/h
d2fv1f1 = (v2xpf1 - (2*v2x0f1) + v2xnf1)/(h**2)
dlf1 = 0.5*(((1/(xd[0]**3))*dfv1f1) - ((1/(xd[0]**2))*d2fv1f1))
dm2f1 = -d2fv1f1 - (3*dlf1*xd[0]**2)

v2x0f2 = uv2x[0]+uv2f2[0]
v2xpf2 = uv2x[1]+uv2f2[1]
v2xnf2 = uv2x[2]+uv2f2[2]
dfv1f2 = (v2xpf2-v2x0f2)/h
d2fv1f2 = (v2xpf2 - (2*v2x0f2) + v2xnf2)/(h**2)
dlf2 = 0.5*(((1/(xd[0]**3))*dfv1f2) - ((1/(xd[0]**2))*d2fv1f2))
dm2f2 = -d2fv1f2 - (3*dlf2*xd[0]**2)

v2x0f3 = uv2x[0]+uv2f3[0]
v2xpf3 = uv2x[1]+uv2f3[1]
v2xnf3 = uv2x[2]+uv2f3[2]
dfv1f3 = (v2xpf3-v2x0f3)/h
d2fv1f3 = (v2xpf3 - (2*v2x0f3) + v2xnf3)/(h**2)
dlf3 = 0.5*(((1/(xd[0]**3))*dfv1f3) - ((1/(xd[0]**2))*d2fv1f3))
dm2f3 = -d2fv1f3 - (3*dlf3*xd[0]**2)

v2x0f4 = uv2x[0]+uv2f4[0]
v2xpf4 = uv2x[1]+uv2f4[1]
v2xnf4 = uv2x[2]+uv2f4[2]
dfv1f4 = (v2xpf4-v2x0f4)/h
d2fv1f4 = (v2xpf4 - (2*v2x0f4) + v2xnf4)/(h**2)
dlf4 = 0.5*(((1/(xd[0]**3))*dfv1f4) - ((1/(xd[0]**2))*d2fv1f4))
dm2f4 = -d2fv1f4 - (3*dlf4*xd[0]**2)

v2x0f5 = uv2x[0]+uv2f5[0]
v2xpf5 = uv2x[1]+uv2f5[1]
v2xnf5 = uv2x[2]+uv2f5[2]
dfv1f5 = (v2xpf5-v2x0f5)/h
d2fv1f5 = (v2xpf5 - (2*v2x0f5) + v2xnf5)/(h**2)
dlf5 = 0.5*(((1/(xd[0]**3))*dfv1f5) - ((1/(xd[0]**2))*d2fv1f5))
dm2f5 = -d2fv1f5 - (3*dlf5*xd[0]**2)

v2f1 = []
v2f2 = []
v2f3 = []
v2f4 = []
v2f5 = []


for i in range(0,j):
    v2f1.append(((dm2f1/2)*x[i]**2) + ((dlf1/4)*x[i]**4))
    v2f2.append(((dm2f2 / 2) * x[i] ** 2) + ((dlf2 / 4) * x[i] ** 4))
    v2f3.append(((dm2f3 / 2) * x[i] ** 2) + ((dlf3 / 4) * x[i] ** 4))
    v2f4.append(((dm2f4 / 2) * x[i] ** 2) + ((dlf4 / 4) * x[i] ** 4))
    v2f5.append(((dm2f5 / 2) * x[i] ** 2) + ((dlf5 / 4) * x[i] ** 4))

vloopf1 = []
vloopf2 = []
vloopf3 = []
vloopf4 = []
vloopf5 = []

for i in range(0,j):
    vloopf1.append(v0[i] + v1x[i] + v1f1[i] + v2f1[i])
    vloopf2.append(v0[i] + v1x[i] + v1f2[i] + v2f2[i])
    vloopf3.append(v0[i] + v1x[i] + v1f3[i] + v2f3[i])
    vloopf4.append(v0[i] + v1x[i] + v1f4[i] + v2f4[i])
    vloopf5.append(v0[i] + v1x[i] + v1f5[i] + v2f5[i])
Tx = 197.871

def ix(r):
    return (r**2) * np.log(1-np.exp(-np.sqrt((r**2)+(x[i]/Tx)**2)))
def if1(r):
    return (r**2) * np.log(1-np.exp(-np.sqrt((r**2)+((mf[1]*x[i])/Tx)**2)))
def if2(r):
    return (r**2) * np.log(1-np.exp(-np.sqrt((r**2)+((mf[2]*x[i])/Tx)**2)))
def if3(r):
    return (r**2) * np.log(1-np.exp(-np.sqrt((r**2)+((mf[3]*x[i])/Tx)**2)))
def if4(r):
    return (r**2) * np.log(1-np.exp(-np.sqrt((r**2)+((mf[4]*x[i])/Tx)**2)))
def if5(r):
    return (r**2) * np.log(1-np.exp(-np.sqrt((r**2)+((mf[5]*x[i])/Tx)**2)))
intx = []
intf1 = []
intf2 = []
intf3 = []
intf4 = []
intf5 = []
for i in range(0,j):
    intx.append(sci.integrate.quad(ix,   0, np.inf))
    intf1.append(sci.integrate.quad(if1, 0, np.inf))
    intf2.append(sci.integrate.quad(if2, 0, np.inf))
    intf3.append(sci.integrate.quad(if3, 0, np.inf))
    intf4.append(sci.integrate.quad(if4, 0, np.inf))
    intf5.append(sci.integrate.quad(if5, 0, np.inf))

uv3x = np.zeros(j)
uv3f1 = np.zeros(j)
uv3f2 = np.zeros(j)
uv3f3 = np.zeros(j)
uv3f4 = np.zeros(j)
uv3f5 = np.zeros(j)
ex = np.zeros(j)
ef1 = np.zeros(j)
ef2 = np.zeros(j)
ef3 = np.zeros(j)
ef4 = np.zeros(j)
ef5 = np.zeros(j)
for i in range(0,j):
    uv3x[i], ex[i] = intx[i]
    uv3f1[i], ef1[i] = intf1[i]
    uv3f2[i], ef2[i] = intf2[i]
    uv3f3[i], ef3[i] = intf3[i]
    uv3f4[i], ef4[i] = intf4[i]
    uv3f5[i], ef5[i] = intf5[i]


v3x = dox*((Tx**4)/(2*np.pi**2))*uv3x
nv3f1 = dof[1]*((Tx**4)/(2*np.pi**2))*uv3f1
nv3f2 = dof[2]*((Tx**4)/(2*np.pi**2))*uv3f2
nv3f3 = dof[3]*((Tx**4)/(2*np.pi**2))*uv3f3
nv3f4 = dof[4]*((Tx**4)/(2*np.pi**2))*uv3f4
nv3f5 = dof[5]*((Tx**4)/(2*np.pi**2))*uv3f5

v3f1 = []
v3f2 = []
v3f3 = []
v3f4 = []
v3f5 = []

for i in range(0,j):
    v3f1.append(v3x[i] + nv3f1[i])
    v3f2.append(v3x[i] + nv3f2[i])
    v3f3.append(v3x[i] + nv3f3[i])
    v3f4.append(v3x[i] + nv3f4[i])
    v3f5.append(v3x[i] + nv3f5[i])

vtotf1 = []
vtotf2 = []
vtotf3 = []
vtotf4 = []
vtotf5 = []

for i in range(0,j):
    vtotf1.append(vloopf1[i] + v3f1[i])
    vtotf2.append(vloopf2[i] + v3f2[i])
    vtotf3.append(vloopf3[i] + v3f3[i])
    vtotf4.append(vloopf4[i] + v3f4[i])
    vtotf5.append(vloopf5[i] + v3f5[i])

#pylab.plot(vloopf1, label='FUCK1')
#pylab.plot(vtotf1 - vtotf1[5], label='FUCKTOT1')
#pylab.legend()
#pylab.xlim(220,250)
#pylab.ylim(-10000,10000)
#pylab.show()


v_t = np.array([((237/v)/197.87), ((330/v)/147.85783), ((360/v)/134.65), ((385/v)/127.87)])

pylab.plot(kf[1:5], v_t)
pylab.show()



