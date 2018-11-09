import numpy as np
import matplotlib.pyplot as plt

filename = "cycleComparison.dat"
dat = np.loadtxt(filename)
anal = np.loadtxt("analytical.dat")


E = dat[:, 1]
M = dat[:, 5]
Cv = dat[:, 7]
chi = dat[:, 6]
cyc = dat[:, -1]

vals = np.array([E, M, Cv, chi])*4.0

print E
print M
print Cv
print chi





plt.show()
