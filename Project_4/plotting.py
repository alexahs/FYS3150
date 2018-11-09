import numpy as np
import matplotlib.pyplot as plt
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


data1 = np.loadtxt("1acceptsVsMC.dat")
data2 = np.loadtxt("2acceptsVsMC.dat")
data3 = np.loadtxt("3acceptsVsMC.dat")
data4 = np.loadtxt("4acceptsVsMC.dat")
data5 = np.loadtxt("5acceptsVsMC.dat")
data6 = np.loadtxt("6acceptsVsMC.dat")
data7 = np.loadtxt("7acceptsVsMC.dat")
data8 = np.loadtxt("8acceptsVsMC.dat")
data9 = np.loadtxt("9acceptsVsMC.dat")
data10 = np.loadtxt("10acceptsVsMC.dat")

allData = np.array([data1, data2, data3, data4, data5, data6, data7, data8, data9, data10])

avgData = np.zeros(99)

for i, j in zip( range(len(data1)), range(len(data1))):
    avgData[i] = sum(allData[j])




print avgData[0]
