## THIS CODE IS IN PYTHON LANGUAGE ##

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# UTILITIES FUNCTION #

def const_fit(x,a):
    const_array = np.linspace(np.min(x), np.max(x), np.size(x))
    for i in range(np.size(x)):
        const_array[i]=a
    return const_array
###

"""
RECO PART
"""

mean_reco = []
sigma_reco = []
for i in range(1,6):
    data_reco = np.loadtxt("10e_"+str(i)+"GeV3_sigma_MR_v2.txt", skiprows=1, delimiter=";")
    mean_reco.append(np.mean(data_reco[:,0]))
    sigma_reco.append(np.std(data_reco[:,0]))

Beam_energy=np.array([1,2,3,4,5])


popt_const_reco, cormat_const_reco = curve_fit(const_fit, Beam_energy,mean_reco, 10)

sigma_const_reco = np.sqrt(cormat_const_reco[0,0])
const_reco = const_fit(Beam_energy, *popt_const_reco)[0]

plt.figure(6)
plt.xlim(0,6)
plt.errorbar(Beam_energy, mean_reco, sigma_reco, marker='o', color='k', linestyle='', label="Reconstruction")
plt.hlines(const_reco, -1, 7, 'r', '-', label="Best constant fit $R_M^{reco}"+f"= {popt_const_reco[0]:.2f}$ cm")
# plt.fill_between([-1,7], const_reco-sigma_const_reco, const_reco+sigma_const_reco, facecolor='g', alpha=0.2)
# plt.fill_between([-1,7], const_reco-2*sigma_const_reco, const_reco+2*sigma_const_reco, facecolor='g', alpha=0.2)
# plt.axhline(const_reco-2*sigma_const_reco, color='g', linestyle='--')
# plt.axhline(const_reco+2*sigma_const_reco, color='g', linestyle='--')

"""
GEANT 4 PART
"""
mean_sim = []
sigma_sim = []
for i in range(1,6):
    data = np.loadtxt("../../sim/better_calculation/10e_"+str(i)+"GeV3_sigma_MR_G4_v2.txt", skiprows=1, delimiter=";")
    mean_sim.append(np.mean(data[:,0]))
    sigma_sim.append(np.std(data[:,0]))

Beam_energy=np.array([1,2,3,4,5])

popt_const_sim, cormat_const_sim = curve_fit(const_fit, Beam_energy,mean_sim, 10)

sigma_const_sim = np.sqrt(cormat_const_sim[0,0])
const_sim = const_fit(Beam_energy, *popt_const_sim)[0]

# plt.figure(7)
plt.xlim(0,6)
plt.errorbar(Beam_energy, mean_sim, sigma_sim, marker='o', color='purple', linestyle='', label="Simulation")
plt.hlines(const_sim, -1, 7, 'coral', '-', label="Best constant fit $R_M^{sim}"+f"= {popt_const_sim[0]:.2f}$ cm")
# plt.hlines(9.04, -1, 7, 'blue', 'dashdot', label="Theoretical value $R_M = 9.04$ cm")
# plt.fill_between([-1,7], const_sim-sigma_const_sim, const_sim+sigma_const_sim, facecolor='y', alpha=0.2)
# plt.fill_between([-1,7], const_sim-2*sigma_const_sim, const_sim+2*sigma_const_sim, facecolor='y', alpha=0.2)
# plt.axhline(const_sim-2*sigma_const_sim, color='y', linestyle='--')
# plt.axhline(const_sim+2*sigma_const_sim, color='y', linestyle='--')
plt.hlines(9.04, -1, 7, 'b', 'dashdot', label="Expected value $R_M = 9.04$ cm")

plt.grid()
plt.legend(loc='best', fontsize=8.5, bbox_to_anchor=(0.3375, 0.5, 0.5, 0.5))
plt.xlabel("Beam Energy [GeV]", fontsize=14)
plt.ylabel("$R_M$ [cm]", fontsize=14)
plt.savefig("Mean_moliere_radius_electrons_rapport.pdf")
plt.savefig("Mean_moliere_radius_electrons_rapport.svg")
