import numpy as np
import scipy as sci
import pylab
from scipy import optimize
from scipy import integrate
from scipy.integrate import odeint

k = 401
lam = 1/8
ms = np.sqrt(lam)
x = np.linspace(0, 1.5, k)
Q = 1
T = 0.61
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

#pylab.plot(x, v0(x))
#pylab.show()

def dU1_dr(u, r):
    return [u[1], -ms**2*u[0] + lam*u[0]**3]


u0 = [-1, 0.0001]
xs1 = np.linspace(-1, 30, k)
Us1 = odeint(dU1_dr, u0, xs1)
ys1 = Us1[:,0]
ysp1 = Us1[:,1]

def fitting(x):
    return(np.tanh(x))



pylab.plot(xs1, fitting((xs1-19)/4))
pylab.plot(xs1, ys1)
pylab.show()

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
    i3b.append(sci.integrate.quad(intb, 0, np.inf))
    i3f.append(sci.integrate.quad(intf, 0, np.inf))
    uv3b[i], eb[i] = i3b[i]
    uv3f[i], ef[i] = i3f[i]
    v3b.append(dofb * ((T**4)/(2*np.pi**2)) * uv3b[i])
    v3f.append(-doff * ((T**4)/(2*np.pi**2)) * uv3f[i])
    v3.append(v3b[i] + v3f[i])

B = vtloop(x[1]) + v3[1]
print(B)

def vtotal(x):
    return vtloop(x) + v3 - B

#pylab.plot(x, vtloop(x), label='vtloop')
#pylab.plot(x, vbloop(x), label='vbloop')
#pylab.plot(x, vfloop(x), label='vfloop')
#pylab.plot(x, v0(x), label='v0')
#pylab.plot(x, vtotal(x), label='Vtotal')
#pylab.legend()
#pylab.show()



#int1 = sci.integrate.quad(vtotal,0, 0.56)
#print(int1)


#Differentiating the LHS


dvtotal = []


for i in range(2, k-1):
    dvtotal.append((vtotal(x)[i+1] - vtotal(x)[i-1])/(2*0.01))


#Fitting to a function

#make in order
def arf(x, a, c, d, f, g):
    return(a*np.sin(x) + c*x**3 + d*x**2 + f*x + g)


arfa, covarfa = sci.optimize.curve_fit(arf, x[1:k-2], dvtotal)

#pylab.plot(dvtotal)
#pylab.plot(arfa[0]*np.sin(x) + arfa[2]*x**3 + arfa[3]*x**2 + arfa[4]*x + arfa[5])
#pylab.show()


#Solving DE


cb = 13.94*16
cf = 13.94
D = -ms**2 + dm2b**2 + (1/12 * T**2 * mx**2 - 1/(6*np.pi)*T*mx) + 1/24 * T**2 * mf**2


def dU_dr(U, r):
    return [U[1], (-2/(r+0.001))*(U[1]) + arfa[0]*np.sin(U[0]) + arfa[1]*U[0]**3 + arfa[2]*U[0]**2 + arfa[3]*U[0] + arfa[4]]


U0 = [0.6074943, 0.0001]
xs = np.linspace(0, 150, k)
Us = odeint(dU_dr, U0, xs)
ys = Us[:,0]
ysp = Us[:,1]



#Integration for S3 graph


def rad(r):
    return np.pi * 4 * r**2


intrad = sci.integrate.quad(rad, 0, 125)

######################################################################################################################


def v0_1(t):
    return (-((ms**2)/2)*t**2) + ((lam/4)*(t**4))


def vbloop_1(t):
    return(((-((ms**2)/2)*t**2) + ((lam/4)*(t**4))) + (A*(dofb*((mx*t)**4))*(np.log(((mx*t)**2)/(Q**2))-1.5)) \
           + (((dm2b / 2) * t ** 2) + ((dlb / 4) * t ** 4)))


def vfloop_1(t):
    return ((-((ms**2)/2)*t**2) + ((lam/4)*(t**4))) + (-A*((doff*((mf*t)**4))*(np.log(((mf*t)**2)/(Q**2))-1.5)))\
           + (((dm2f / 2) * t ** 2) + ((dlf / 4) * t ** 4))


def vtloop_1(t):
    return vbloop_1(t) + vfloop_1(t) - v0_1(t)


def intb_1(r):
    return (r**2) * np.log(1-np.exp(-np.sqrt((r**2) + ((mx*ys[i])/T)**2)))


def intf_1(r):
    return (r**2) * np.log(1+np.exp(-np.sqrt((r**2) + ((mf*ys[i])/T)**2)))


i3b_1 = []
i3f_1 = []
uv3b_1 = np.zeros(k)
uv3f_1 = np.zeros(k)
eb_1 = np.zeros(k)
ef_1 = np.zeros(k)
v3b_1 = []
v3f_1 = []
v3_1 = []
for i in range(0, k):
    i3b_1.append(sci.integrate.quad(intb_1,0,np.inf))
    i3f_1.append(sci.integrate.quad(intf_1,0,np.inf))
    uv3b_1[i], eb_1[i] = i3b_1[i]
    uv3f_1[i], ef_1[i] = i3f_1[i]
    v3b_1.append(dofb * ((T**4)/(2*np.pi**2)) * uv3b_1[i])
    v3f_1.append(-doff * ((T**4)/(2*np.pi**2)) * uv3f_1[i])
    v3_1.append(v3b_1[i] + v3f_1[i])

B_1 = vtloop_1(ys[1]) + v3_1[1]


def vtotal_1(ys):
    return vtloop_1(ys) + v3_1 - B_1

def integrandofr(r):
    return r**2


integrand = []
integrand.append((vtotal_1(ys) + 0.5*ysp**2)*integrandofr(xs))


ints3 = []
ints3 = np.trapz(integrand, xs)

print(4*np.pi*ints3/T)


def hypt(xs):
    return (U0[0]/2)*(1-np.tanh((xs-80)/23))


#print(vtotal(q[0]))

#pylab.plot(vtotal_1(ys))


#print(ys)
#pylab.plot(xs, ys)
#pylab.xlabel('r')
#pylab.ylabel('Phi')
#pylab.plot(xs, hypt(xs))
#pylab.show()

#pylab.plot(xs, ysp**2 / 2)
#pylab.plot(xs, ((vtotal_1(ys)) + ysp**2 / 2)*integrandofr(xs))
#pylab.show()

#print(vtotal_1(ys)[0])
#print(ys[0])

#pylab.plot(xs[1:400], 2* ysp[1:400]**2/xs[1:400])
#pylab.show()

#p = 2 * ysp[1:400]**2/xs[1:400]
#xs2 = np.linspace(1, 150, k-2)
#pp = sci.integrate.trapz(p, xs2, k-2)
#print(pp)
#print(vtotal_1(ys)[0])
#print(vtotal_1(ys)[1]/pp)