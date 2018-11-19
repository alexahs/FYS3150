import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab


from scipy.stats import norm
from scipy.integrate import quad
plt.style.use('ggplot')
# import seaborn as sns
# sns.set()
# sns.set_style('white')
# plt.xkcd();
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')


# ===============OPPGAVE 4B========================
# filename = "cycleComparison.txt"
# dat = np.loadtxt(filename)
# anal = np.loadtxt("analytical.txt")
#
#
# E = dat[:, 1]
# M = dat[:, 5]
# Cv = dat[:, 7]
# chi = dat[:, 6]
# cyc = dat[:, -1]
#
# anal = anal/4.0
#
# # errors = np.zeros((len(E), 4))
# errE = np.zeros(len(E))
# errM = np.zeros(len(E))
# errCv = np.zeros(len(E))
# errChi = np.zeros(len(E))
#
#
#
# # vals = np.array([E, M, Cv, chi])*4.0
#
# def absE(A, S):
#     # print A
#     # print S
#     return abs((A - S)/A)
#
# errE = absE(anal[0], E)
# errM = absE(anal[4], M)
# errCv = absE(anal[5], Cv)
# errChi = absE(anal[6], chi)
#
# figname = "errors_E_M.pdf"
# # plt.loglog(cyc, errE)
# # plt.loglog(cyc, errM)
# plt.loglog(cyc, errE)
# plt.loglog(cyc, errM)
# plt.legend([r"$\langle E \rangle$", r"$|M|$"])
# plt.xlabel("No. of MC Cycles")
# plt.ylabel("Absolute error")
# plt.savefig(figname)
#
#
# plt.show()
# ===================================================

# # =============OPPGAVE 4C============================
# cyc = 200
# cyc_string = str(cyc)
#
# ER_filename = "0calibrate20Cycles"+ cyc_string +"Ordered0.bin"
# EO_filename = "0calibrate20Cycles"+ cyc_string +"Ordered1.bin"
# MR_filename = "1calibrate20Cycles"+ cyc_string +"Ordered0.bin"
# MO_filename = "1calibrate20Cycles"+ cyc_string +"Ordered1.bin"
#
# cycles = np.linspace(1, cyc, cyc+2)
# ERvalues = np.fromfile(ER_filename, dtype=np.int32)/400.0
# EOvalues = np.fromfile(EO_filename, dtype=np.int32)/400.0
# MRvalues = np.fromfile(MR_filename, dtype=np.int32)/400.0
# MOvalues = np.fromfile(MO_filename, dtype=np.int32)/400.0
#
#
# meanE = (sum(ERvalues) + sum(EOvalues))/float(len(ERvalues) + len(EOvalues))
# meanM = (sum(MRvalues) + sum(MOvalues))/float(len(MRvalues) + len(MOvalues))
#
# print meanE, meanM
#
# plt.plot(cycles, ERvalues)
# plt.plot(cycles, EOvalues)
# plt.plot((1, cyc), (meanE, meanE), "k")
# plt.xlabel("No. of MC cycles")
# plt.ylabel(r"$\langle E \rangle$")
# plt.legend(["Random lattice", "Ordered lattice", r"$\overline{E}$"])
# # plt.axis([0, cyc+2, -2.1, -1.5 ])
# plt.savefig("equilibriumTime_E_T24.pdf")
# plt.show()
# #######################
# plt.plot(cycles, MRvalues)
# plt.plot(cycles, MOvalues)
# plt.plot((1, cyc), (meanM, meanM), "k")
# plt.xlabel("No. of MC cycles")
# plt.ylabel(r"$\langle |M| \rangle$")
# plt.legend(["Random lattice", "Ordered lattice", r"$\overline{|M|}$"])
# # plt.axis([0, cyc+2, -2.1, -1.5 ])
# plt.savefig("equilibriumTime_M_T24.pdf")
# plt.show()
# ===================================================


# dataOrder = np.loadtxt("acceptsRatioVsT.txt")
# dataRandom = np.loadtxt("acceptsRatioVsTdisordered.txt")
#
# T = dataOrder[:, 0]
# ratOrder = dataOrder[:, 1]
# ratRandom = dataRandom[:, 1]
#
# plt.plot(T, ratOrder)
# plt.plot(T, ratRandom)
#
# plt.xlabel("Temperature [kT/J]")
# plt.ylabel("Ratio of accepted configurations")
# plt.legend(["Ordered lattice", "Random lattice"])
# plt.savefig("acceptRatio.pdf")
# plt.show()

# ====================OPPGAVE 4E=================================

# data40 = np.loadtxt("dim40CriticalTemps100000Cycles.dat", skiprows = 1)
# data60 = np.loadtxt("dim60CriticalTemps100000Cycles.dat", skiprows = 1)
# data80 = np.loadtxt("dim80CriticalTemps100000Cycles.dat", skiprows = 1)
# data100 = np.loadtxt("dim100CriticalTemps100000Cycles.dat", skiprows = 1)
#
# temps = data40[:,0]
# data40 = data40[:,1:]
# data60 = data60[:,1:]
# data80 = data80[:,1:]
# data100 = data100[:,1:]
#
#
# plt.plot(temps, data40[:,3])
# plt.plot(temps, data60[:,3])
# plt.plot(temps, data80[:,3])
# plt.plot(temps, data100[:,3])
# plt.xlabel("Temperature [kT/J]")
# plt.ylabel(r"$C_v$")
# plt.legend("40 x 40", "60 x 60", "80 x 80", "100 x 100")
# plt.savefig("Heatcap.pdf")
# plt.show()
#
# plt.plot(temps, data40[:,2])
# plt.plot(temps, data60[:,2])
# plt.plot(temps, data80[:,2])
# plt.plot(temps, data100[:,2])
# plt.xlabel("Temperature [kT/J]")
# plt.ylabel(r"$\chi$")
# plt.savefig("Susceptibility.pdf")
# plt.show()
#
# plt.plot(temps, data40[:,1])
# plt.plot(temps, data60[:,1])
# plt.plot(temps, data80[:,1])
# plt.plot(temps, data100[:,1])
# plt.xlabel("Temperature [kT/J]")
# plt.ylabel(r"$\langle \overline{M} \rangle$")
# plt.savefig("Magnetization.pdf")
# plt.show()
#
# plt.plot(temps, data40[:,0])
# plt.plot(temps, data60[:,0])
# plt.plot(temps, data80[:,0])
# plt.plot(temps, data100[:,0])
# plt.ylabel(r"$\langle E \rangle$ ")
# plt.xlabel("Temperature [kT/J]")
# plt.savefig("Energy.pdf")
# plt.show()



# =================OPPGAVE 4D==============================
# E = np.fromfile("energyProbabilityT2.400000Cycles100000.bin", dtype=np.int32)
# print E
# indices = int(round(len(E)*0.9))
# E = E[indices:]
# a = -800
# b = -300
#
# x = np.linspace(a, b, 1000)
# num_bins = 50
# mu, sigma = norm.fit(E)
#
# n, bins, patches = plt.hist(E, num_bins, density=True, facecolor='blue', align='mid', alpha=0.75)
# y = mlab.normpdf( bins, mu, sigma)
#
#
#
# # plt.plot(bins, y, 'r--')
# # E = E/(20*1e5)
# # plt.hist(E, num_bins, density=True, color='b')
# plt.plot(bins, y, color='r')
# plt.xlabel("Lattice energy")
# plt.ylabel("Probability density")
# plt.legend(["Fitted normal distribution", "Sampeled energy"])
# plt.savefig("EnergyProbability.pdf")
# # print min(E[indices:])
# plt.show()

# =========HIGH RES PLOTS===============


data80 = np.loadtxt("HIGHRESdim120CriticalTemps100000Cycles.dat", skiprows = 1)
data100 = np.loadtxt("HIGHRESdim140CriticalTemps100000Cycles.dat", skiprows = 1)

temps = data80[:,0]
data80 = data80[:,1:]
data100 = data100[:,1:]


plt.plot(temps, data80[:,3])
plt.plot(temps, data100[:,3])
plt.xlabel("Temperature [kT/J]")
plt.ylabel(r"$C_v$")
plt.legend(["120 x 120", "140 x 140"])
plt.savefig("Heatcap.pdf")
plt.show()

plt.plot(temps, data80[:,2])
plt.plot(temps, data100[:,2])
plt.xlabel("Temperature [kT/J]")
plt.ylabel(r"$\chi$")
plt.savefig("Susceptibility.pdf")
plt.show()

plt.plot(temps, data80[:,1])
plt.plot(temps, data100[:,1])
plt.xlabel("Temperature [kT/J]")
plt.ylabel(r"$\langle \overline{M} \rangle$")
plt.savefig("Magnetization.pdf")
plt.show()

plt.plot(temps, data80[:,0])
plt.plot(temps, data100[:,0])
plt.ylabel(r"$\langle E \rangle$ ")
plt.xlabel("Temperature [kT/J]")
plt.savefig("Energy.pdf")
plt.show()
