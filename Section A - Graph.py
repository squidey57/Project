import numpy as np
import pylab as pl
import scipy as sci

#For Nx=1 the system works best as there are not too many particles in the system with too high masses

kx = [1.00, 1.25, 1.50, 1.75, 1.90]
Tc = [0.6175, 0.548, 0.48, 0.39, 0.275]
Vc = [112/200, 149/200, 176/200, 192/200, 199/200]
R = [0.907, 1.369, 1.833, 2.46, 3.618]

pl.plot(kx, R)
pl.xlabel('Kx')
pl.ylabel('Vc/Tc')
pl.show()