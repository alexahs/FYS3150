import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')



# fname = "results_r.txt"
#
# data = np.loadtxt(fname)
#
# rho = data[:, 0]
# psi = abs(data[:, 23])**2
#
# print abs(data[:, 23])
# print data[:, 1]
#
# plt.plot(rho, psi)


###JACOBI ITERATIONS###
# n = np.array([10, 25, 50, 100, 150, 200, 300, 400])
# iterations = np.array([135,1110, 4280, 17600, 39900, 71100, 160700, 287100])
# x = np.linspace(10, 400, 400)
# y = lambda t : t**2
# plt.style.use('ggplot')
# plt.semilogy(x, y(x))
# plt.semilogy(n, iterations)
# plt.xlabel("Number of grid points")
# plt.ylabel("Number of iterations")
# plt.legend(["$n^2$", "Jacobi iterations"])

###RUN TIMES###
# n = np.array([10, 25, 50, 100, 150, 200, 300, 400])
# times = [2.4e-4, 4.5e-3, 2.4e-2, 0.35, 1.67, 5.153, 25.4, 80]
# plt.style.use('ggplot')
# plt.plot(n, times)
# plt.xlabel("Number of grid points")
# plt.ylabel("Number of iterations")
# plt.legend(["$n^2$", "Jacobi iterations"])



###ANALYTICAL BUCKLING BEAM###

# 
# data = np.loadtxt("results_b.txt")
# rho = data[:, 0]
# comp = data[:, 1]
# anal = data[:, 2]
#
# plt.plot(rho, comp)
# plt.plot(rho, anal)




###RHOMAX###
# rhoMax = np.linspace(2, 6, 17)
# ratio0 = np.array([0.857274, 0.929333, 0.969761, 0.988908, 0.996536, 0.999092, 0.999817, 0.999994,
# 1.00003, 1.00005, 1.00005, 1.00006, 1.00007, 1.00007, 1.00008, 1.00009, 1.00009])
# ratio0 = 1.0/ratio0
# ratio1 = [1.57141,1.33307,1.18358,1.0932,1.04236,1.01674,1.00559,1.00152,1.00029,0.999974,0.999899,0.999875,0.999859,0.999845,0.999829,0.999814,0.999797,
# ]
# ratio2 = [2.10123,1.7195,1.46064,1.28397,1.16506,1.08815,1.04194,1.01717,1.00583,1.00154,1.00022,0.999879,0.999793,0.99976,0.999735,0.999711,0.999685,
# ]
# ratio3 = [2.66932,2.1525,1.79302,1.53768,1.35467,1.22402,1.13272,1.07178,1.03424,1.01383,1.00448,1.00102,0.999986,0.999723,0.999648,0.999608,0.999572]
# plt.style.use('ggplot')
# plt.plot(rhoMax, ratio0)
# plt.plot(rhoMax, ratio1)
# plt.plot(rhoMax, ratio2)
# plt.plot(rhoMax, ratio3)
# plt.legend(["$\lambda = 3$", "$\lambda = 7$", "$\lambda = 11$", "$\lambda = 15$"])
# plt.xlabel(r"$\rho_{max}$")
# plt.ylabel(r"$\lambda /\lambda_{computed}$")

##N VALUES###
# nvals = [10, 30, 50, 70, 90, 110, 130, 150, 170, 190, 210, 230, 250, 270, 290, 310, 330, 350, 370, 385, 400]
# data = np.loadtxt("analyze_n.txt")
#
#
# plt.style.use('ggplot')
# plt.plot(nvals, data[:, 0])
# plt.plot(nvals, data[:, 1])
# plt.plot(nvals, data[:, 2])
# plt.plot(nvals, data[:, 3])
# plt.xlabel("No. of grid points")
# plt.ylabel(r"$\lambda /\lambda_{computed}$")
# plt.legend(["$\lambda = 3$", "$\lambda = 7$", "$\lambda = 11$", "$\lambda = 15$"], loc='best')


###OMEGA R###
omega_r001 = np.loadtxt("results_omega_0_01.txt")
omega_r05 = np.loadtxt("results_omega_0_5.txt")
omega_r1 = np.loadtxt("results_omega_1.txt")
omega_r5 = np.loadtxt("results_omega_5.txt")

plt.style.use('ggplot')
plt.plot(omega_r001[:, 0], omega_r001[:, 1]**2)
plt.plot(omega_r05[:, 0], omega_r05[:, 1]**2)
plt.plot(omega_r1[:, 0], omega_r1[:, 1]**2)
plt.plot(omega_r5[:, 0], omega_r5[:, 1]**2)
plt.legend(["$\omega_r = 0.01, \lambda_0 = 0.828$", "$\omega_r = 0.5,  \lambda_0 = 2.230$", "$\omega_r = 1,   \lambda_0 = 4.058$", "$\omega_r = 5,  \lambda_0 = 17.444$"], loc='best')
plt.xlabel(r"$\rho$")
plt.ylabel(r"Probability density $|\psi|^2$")






plt.show()
