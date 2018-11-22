import numpy as np

x =np.linspace(0,10,11)

def v0(x):
    return x**2 + 5

A = np.array([1,1,1,1,1,1,1,1,1,1,1])

print(v0(x) + A)
