import numpy as np
import scipy as sci
import pylab
from scipy import optimize
from scipy import integrate

l = 1 / 8
m = np.sqrt(1 / 8)
k = 1200
x = np.linspace(0, 6, k)
t = np.ones(k)
mx = 10
dox = 3
A = (1 / (64 * np.pi ** 2))
Q = 1
dof = 4
h = 0.01
xd = np.array([1, (1 + h), (1 - h)])

uv2b = []
uv2f = []

for i in range(0, 3):
    uv2b.append(A * dox * ((3 * xd[i]) ** 4) * (np.log(((3 * xd[i]) ** 2) / Q ** 2) - 1.5))
    uv2f.append(-A * dof * ((1 * xd[i]) ** 4) * (np.log(((1 * xd[i]) ** 2) / Q ** 2) - 1.5))

v2x0b1 = uv2b[0]
v2xpb1 = uv2b[1]
v2xnb1 = uv2b[2]
dfv1b1 = ((v2xpb1 - v2xnb1) / (2 * h))
d2fv1b1 = ((v2xpb1 - (2 * v2x0b1) + v2xnb1) / (h ** 2))
dlb = 0.5 * (dfv1b1 - d2fv1b1)
dm2b = 0.5 * (1 / xd[0]) * ((d2fv1b1 * xd[0]) - (3 * dfv1b1))

print(d2fv1b1 + dm2b + (3*dlb) )
v2x0f1 = uv2f[0]
v2xpf1 = uv2f[1]
v2xnf1 = uv2f[2]
dfv1f1 = (v2xpf1 - v2xnf1) / (2 * h)
d2fv1f1 = (v2xpf1 - (2 * v2x0f1) + v2xnf1) / (h ** 2)
dlf = 0.5 * (dfv1f1 - d2fv1f1)
dm2f = 0.5 * (1 / xd[0]) * ((d2fv1f1 * xd[0]) - (3 * dfv1f1))

def vnloop(x):
    return (-((m ** 2) / 2) * x ** 2) + ((l / 4) * x ** 4) + \
           (A * dox * ((3 * x) ** 4) * (np.log(((3 * x ** 2) / Q ** 2)) - 1.5)) + \
           (-A * dof * ((1 * x) ** 4) * (np.log(((1 * x ** 2) / Q ** 2)) - 1.5)) + \
           (((dm2b / 2) * x ** 2) + ((dlb / 4) * x ** 4)) + (((dm2f / 2) * x ** 2) + ((dlf / 4) * x ** 4))

T = 2.0078
def intb(r):
    return (r**2) * np.log(1-np.exp(-np.sqrt((r**2) + ((3*x[i])/T)**2)))

def intf(r):
    return (r**2) * np.log(1+np.exp(-np.sqrt((r**2) + ((1*x[i])/T)**2)))

i3b = []
i3f = []
uv3b = np.zeros(k)
uv3f = np.zeros(k)
eb = np.zeros(k)
ef = np.zeros(k)
v3b = []
v3f = []
v3 = []
for i in range(0,k):
    i3b.append(sci.integrate.quad(intb,0,np.inf))
    i3f.append(sci.integrate.quad(intf,0,np.inf))
    uv3b[i], eb[i] = i3b[i]
    uv3f[i], ef[i] = i3f[i]
    v3b.append(dox * ((T**4)/(2*np.pi**2)) * uv3b[i])
    v3f.append(-dof * ((T**4)/(2*np.pi**2)) * uv3f[i])
    v3.append(v3b[i] + v3f[i])

Vtot = vnloop(x) + v3
Vtot1 = vnloop(x[1]) + v3[1]
print(np.min(Vtot - Vtot1))
pylab.plot(vnloop(x) , label = 'VnLoop')
pylab.plot(Vtot - Vtot1, label = 'VTotal')
pylab.ylim(-0.001,0.001)
pylab.xlim(475,500)
pylab.legend()
pylab.show()


v_tf = np.array([(112.5/200)/(0.6175), (48/200)/(0.7572), (20/200)/(1.1125), (10/200)/1.5035, (8/200)/1.8994])
v_tb = np.array([(113/200)/(0.6175), (337/200)/0.815, (490/200)/2.008, (569/200)/3.782, (672/200)/5.02])
k = np.array([1,2,3,4,5])

#pylab.plot(k, v_tf, label = 'Fermion Mass Variation')
#pylab.plot(k, v_tb, label = 'Boson Mass Variation')
#pylab.legend()
#pylab.savefig('boobaa')
#pylab.show()


