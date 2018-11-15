import numpy as np
import scipy as sci
import pylab
from scipy import integrate
from scipy import optimize


l = 1/8
j = 1000
v = 256
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

pylab.plot(vloopf1, label='FUCK1')
pylab.plot(vloopf2, label='FUCK2')
pylab.plot(vloopf3, label='FUCK3')
pylab.plot(vloopf4, label='FUCK4')
pylab.plot(vloopf5, label='FUCK5')
pylab.legend()
pylab.show()



