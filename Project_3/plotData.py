import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.signal import find_peaks
plt.style.use('ggplot')
convert = 3600*180/np.pi
#
data = np.loadtxt("verlet.dat")
# ###3 DIMENSIONS###
n = len(data[:,0])
p_no = int(len(data[1,:])/3.0)

planetlist = ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]


w = 13; h = 8
fig = plt.figure(figsize=(w,h))
ax = fig.add_subplot(111, projection='3d')

for i in range(p_no):
    plt.plot(data[:, i*3], data[:, i*3+1], data[:,i*3+2])

ax.set_xlabel("x [AU]")
ax.set_ylabel("y [AU]")
ax.set_zlabel("z [AU]")
plt.legend(planetlist)
plt.show()


# plt.plot(data[:,0], data[:,1], data[:,2])
# plt.plot(data[:,3], data[:,4], data[:,5])
# plt.plot(data[:,6], data[:,7], data[:,8])
# plt.plot(data[:,9], data[:,10], data[:,11])
# plt.plot(data[:,12], data[:,13], data[:,14])
# plt.plot(data[:,15], data[:,16], data[:,17])
# plt.plot(data[:,18], data[:,19], data[:,20])


# plt.show()

###2 DIMENSIONS###

# n = len(data[:,0])
# p_no = int(len(data[1,:])/3.0)
#
# # for i in range(p_no):
# #     plt.plot(data[:, i*3], data[:, i*3+1])
#
# # plt.plot(data[:, 0], data[:, 1], 'b-')
# plt.plot(data[:, 3], data[:, 4], 'r-')
# # plt.plot(data[:, 6], data[:, 7], 'r-')
#
#
# # plt.plot(data[-1, 0], data[-1, 1], 'bo')
#
# plt.plot(0,0, 'yo')
# plt.plot(data[-1, 3], data[-1, 4], 'bo')
# # plt.plot(data[-1, 6], data[-1, 7], 'ro')
# plt.xlabel("x [AU]")
# plt.ylabel("y [AU]")
# plt.show()

# plt.plot(data[:,0], data[:,1])

# peri_ind, _ = find_peaks(data[:,0])
# npoints =  len(peri_ind)
#
# x = data[:, 1]
# y = data[:, 2]
# t = np.linspace(0, 100, npoints)
#
# print len(data[:, 0])
# peri_angles = np.zeros(npoints)
# for i, j in zip(peri_ind, range(npoints)):
#     peri_angles[j] = np.arctan(y[i]/x[i])*arcsec_convert
#
# plt.plot(t, peri_angles)
# plt.show()

###PRECESSION###
#
# data = np.loadtxt("precession_GR.dat")
#
# alength = len(data[:,0])
# x = data[:,0]
# y = data[:,1]
# t = np.linspace(0, 100, alength)
#
# arcsecs = np.zeros(alength)
# for i in range(alength):
#     arcsecs[i] = np.arctan(y[i]/x[i])*convert
#
#
# p = np.polyfit(t, arcsecs, 1)
#
# b = p[0]
# m = p[1]
# yline = np.polyval(p, t)
# print yline[-1] - yline[0]
#
#
# plt.plot(t, arcsecs)
# plt.plot(t, yline)
# plt.xlabel("Time [years]")
# plt.ylabel("Seconds of arc")
# plt.legend(["Precession", "Fitted line"])
# plt.show()


###2D
#
# data = np.loadtxt("sunearthbeta.dat")
#
# plt.plot(data[:,2], data[:,3])
# plt.plot(0, 0, 'yo')
# plt.show()




###RUN TIMES

# h = np.linspace(1, 8, 8)*-1
# euler = np.array([1.2e-6, 7.6e-6, 6.5e-5, 1.1e-3, 1.0e-2, 6.3e-2, 0.53, 5.2])
# verlet = np.array([3.9e-6, 1.7e-5, 9.5e-5, 1.6e-3, 1.1e-2, 6.5e-2, 0.58, 5.8])
# for i in range(len(h)):
#     h[i] = 10**h[i]
#
# plt.loglog(h, euler)
# plt.loglog(h, verlet)
# plt.xlabel("Step size")
# plt.ylabel("Run time [s]")
# plt.legend(["Euler", "Velocity Verlet"])
# plt.show()
