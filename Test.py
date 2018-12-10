import numpy as np
import pylab as pl


sigma = 4338732.6
dv = 0.02


def s(r):
    return -(np.pi * (4/3) * r**3 * dv) + (4 * np.pi * sigma * r**2)


r = np.linspace(0, 800000000)
pl.plot(s(r))
pl.show()
