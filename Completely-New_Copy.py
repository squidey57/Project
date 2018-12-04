import numpy as np
import scipy as sci
from scipy import integrate
from scipy import optimize
import pylab

#dhsjahjshfhdfhsdjhgdshgjdshghdfhdihfshfdsfd

k = 400
x = np.linspace(0,k,k)
m = 125/np.sqrt(2)

l = 0.129
j = 2
gf = np.array([-12,6,3,1])
m0 = np.array([173,80,91])
gc = np.array([((np.sqrt(2)*m0[0])/246),(2*m0[1]/246), np.sqrt(((4*m0[2]**2)/(246**2))-(2*m0[1]/246)**2)])
m_x = np.array([(gc[0])/np.sqrt(2), (gc[1]/2), np.sqrt((((gc[1]**2)+(gc[2]**2))/4)), j])
Q = 173
A = (1/(64*np.pi**2))

def v0(x):
    return -(((m**2)/2)*x**2)+((l/4)*x**4)


def v1t(x):
    return A*gf[0]*(((m_x[0]*x)**4)*(np.log(((m_x[0]*x)**2)/(Q**2))-1.5))
def v1w(x):
    return A*gf[1]*(((m_x[1]*x)**4)*(np.log(((m_x[1]*x)**2)/(Q**2))-1.5))
def v1z(x):
    return A*gf[2]*(((m_x[2]*x)**4)*(np.log(((m_x[2]*x)**2)/(Q**2))-1.5))
def v1x(x):
    return A*gf[3]*(((m_x[3]*x)**4)*(np.log(((m_x[3]*x)**2)/(Q**2))-1.5))

#Requires the calculation of a first and second derivative
h = 0.001
xd = np.array([246,246+h,246-h])

#It is possible to do the derivatives of V2 algebraically, but V1 is not as easy to calculate
#The masses will be slightly different for +h and -h for each case, so they must calculated seperately

mh0 = np.array([(gc[0]/np.sqrt(2))*xd[0], (gc[1]/2)*xd[0], np.sqrt((((gc[1]**2)+(gc[2]**2))/4)*xd[0]**2), (j*xd[0])])
mhp = np.array([(gc[0]/np.sqrt(2))*xd[1], (gc[1]/2)*xd[1], np.sqrt((((gc[1]**2)+(gc[2]**2))/4)*xd[1]**2), (j*xd[1])])
mhn = np.array([(gc[0]/np.sqrt(2))*xd[2], (gc[1]/2)*xd[2], np.sqrt((((gc[1]**2)+(gc[2]**2))/4)*xd[2]**2), (j*xd[2])])


#Now to calculate a +h and -h for the V1

v10 = []
v1p = []
v1n = []
for i in range(0,4):
    v10.append(gf[i]*(mh0[i]**4)*(np.log((mh0[i]**2)/Q**2)-1.5))
    v1p.append(gf[i] * (mhp[i] ** 4) * (np.log((mhp[i] ** 2) / Q ** 2) - 1.5))
    v1n.append(gf[i] * (mhn[i] ** 4) * (np.log((mhn[i] ** 2) / Q ** 2) - 1.5))

fh0 = (1/(64*np.pi**2))*np.sum(v10)
fhp = (1/(64*np.pi**2))*np.sum(v1p)
fhn = (1/(64*np.pi**2))*np.sum(v1n)

#Now it is possible to the derivates of V1

dfv1 = (fhp-fh0)/h
d2fv1 = (fhp - (2*fh0) + fhn)/(h**2)

#Now that we have these differential values, it is now possible to calculate dm2 and dl

dl = 0.5*(((1/(xd[0]**3))*dfv1) - ((1/(xd[0]**2))*d2fv1))
dm2 = -d2fv1 - (3*dl*xd[0]**2)

def v2(x):
    return ((dm2/2)*x**2) + ((dl/4)*x**4)


def vnloop(x):
    return v0(x) + v1t(x) + v1w(x) + v1z(x) + v2(x) + v1x(x)

T = 105


def it(r):
    return (r**2)*np.log(1+np.exp(-np.sqrt((r**2)+((((m_x[0])*(x[i]))**2)/T**2))))
def iw(r):
    return (r**2)*np.log(1-np.exp(-np.sqrt((r**2)+((((m_x[1])*(x[i])))**2)/T**2)))
def iz(r):
    return (r**2)*np.log(1-np.exp(-np.sqrt((r**2)+((((m_x[2])*(x[i])))**2)/T**2)))
def ix(r):
    return (r**2)*np.log(1-np.exp(-np.sqrt((r**2)+((((m_x[3])*(x[i])))**2)/T**2)))


#Now to do the integral for all possible values of the masses



#From this integral function, each set is given 2 values, one is the actual result, the other is
#the associated error term, it is necessary to seperate these 2 values to continue
intt = []
intw = []
intz = []
intx = []
for i in range(0,k):
    intt.append(sci.integrate.quad(it,0,np.inf))
    intw.append(sci.integrate.quad(iw,0,np.inf))
    intz.append(sci.integrate.quad(iz,0,np.inf))
    intx.append(sci.integrate.quad(ix,0,np.inf))

uv3t = np.zeros(k)
uv3w = np.zeros(k)
uv3z = np.zeros(k)
uv3x = np.zeros(k)
et = np.zeros(k)
ew = np.zeros(k)
ez = np.zeros(k)
ex = np.zeros(k)
for i in range(0,k):
    uv3t[i], et[i] = intt[i]
    uv3w[i], ew[i] = intw[i]
    uv3z[i], ez[i] = intz[i]
    uv3x[i], ex[i] = intx[i]

v3t = gf[0]*((T**4)/(2*np.pi**2))*uv3t
v3w = gf[1]*((T**4)/(2*np.pi**2))*uv3w
v3z = gf[2]*((T**4)/(2*np.pi**2))*uv3z
v3x = gf[3]*((T**4)/(2*np.pi**2))*uv3x

#Now to collate these terms to give a full v3

v3 = []
for i in range(0,k):
    v3.append(v3t[i]+v3w[i]+v3z[i]+v3x[i])

vtot = vnloop(x)+v3
print(v1t(x))
#pylab.plot(vnloop(x), label='ggyugu')
pylab.plot(vtot - vtot[1], label='j=2')
pylab.legend()
#pylab.savefig('Thermal and X boson, j=-4')
#pylab.ylim(-300,1000)
pylab.xlim(0,250)
pylab.show()
