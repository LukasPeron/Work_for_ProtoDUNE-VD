## THIS CODE IS IN PYTHON LANGUAGE ##

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

mean = []
sigma = []
for i in range(1,6):
    data = np.loadtxt("10e_"+str(i)+"GeV3_sigma_RL_G4_v2.txt", skiprows=1, delimiter=";")
    mean.append(np.mean(data[:,0]))
    sigma.append(np.std(data[:,0]))

Beam_energy=np.array([1,2,3,4,5])

def const_fit(x,a):
    const_array = np.linspace(np.min(x), np.max(x), np.size(x))
    for i in range(np.size(x)):
        const_array[i]=a
    return const_array

popt, cormat = curve_fit(const_fit, Beam_energy,mean,1)
popt_const, cormat_const = curve_fit(const_fit, Beam_energy,mean, 10)

sigma_const = np.sqrt(cormat_const[0,0])
const = const_fit(Beam_energy, *popt_const)[0]

plt.figure(6)
plt.xlim(0,6)
plt.errorbar(Beam_energy, mean, sigma, marker='o', color='k', linestyle='', label="Data")
# plt.plot(Beam_energy, fit_fct(Beam_energy, *popt),'-r', label="best affine fit")
# plt.plot(Beam_energy, const_fit(Beam_energy, *popt_const),'-r', label=f"best constant fit $R_M = {popt_const[0]:.2e}$ cm")
plt.hlines(const, -1, 7, 'red', '-', label=f"Best constant fit $X_0 = {popt_const[0]:.2f}$ cm")
plt.hlines(14, -1, 7, 'blue', 'dashdot', label="Known value $X_0 = 14$ cm")

plt.fill_between([-1,7], const-sigma_const, const+sigma_const, facecolor='g', alpha=0.2)
plt.fill_between([-1,7], const-2*sigma_const, const+2*sigma_const, facecolor='g', alpha=0.2)
plt.axhline(const-2*sigma_const, color='g', linestyle='--')
plt.axhline(const+2*sigma_const, color='g', linestyle='--')
plt.grid()
plt.legend(loc='best', bbox_to_anchor=(0.5, -0.15, 0.5, 0.5))
plt.xlabel("Beam Energy [GeV]", fontsize=14)
plt.ylabel("Mean Radiation length for 10 electrons [cm]", fontsize=14)
plt.title("Mean Radiation length for electromagnetic showers\ncause by an electrons for various beam energy\nfrom simulated data with 95% CL", fontsize=14)
plt.savefig("Mean_radiation_length_electrons_G4_v2.svg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

mean = []
sigma = []
for i in range(1,6):
    data = np.loadtxt("10e_"+str(i)+"GeV3_sigma_RL_G4_v2.txt", skiprows=1, delimiter=";")
    mean.append(np.mean(data[:,0]))
    sigma.append(np.std(data[:,0]))

Beam_energy=np.array([1,2,3,4,5])

def const_fit(x,a):
    const_array = np.linspace(np.min(x), np.max(x), np.size(x))
    for i in range(np.size(x)):
        const_array[i]=a
    return const_array

popt, cormat = curve_fit(const_fit, Beam_energy,mean,1)
popt_const, cormat_const = curve_fit(const_fit, Beam_energy,mean, 10)

sigma_const = np.sqrt(cormat_const[0,0])
const = const_fit(Beam_energy, *popt_const)[0]

plt.figure(6)
plt.xlim(0,6)
plt.errorbar(Beam_energy, mean, sigma, marker='o', color='k', linestyle='', label="Data")
# plt.plot(Beam_energy, fit_fct(Beam_energy, *popt),'-r', label="best affine fit")
# plt.plot(Beam_energy, const_fit(Beam_energy, *popt_const),'-r', label=f"best constant fit $R_M = {popt_const[0]:.2e}$ cm")
plt.hlines(const, -1, 7, 'red', '-', label=f"Best constant fit $X_0 = {popt_const[0]:.2f}$ cm")
plt.hlines(14, -1, 7, 'blue', 'dashdot', label="Known value $X_0 = 14$ cm")

plt.fill_between([-1,7], const-sigma_const, const+sigma_const, facecolor='g', alpha=0.2)
plt.fill_between([-1,7], const-2*sigma_const, const+2*sigma_const, facecolor='g', alpha=0.2)
plt.axhline(const-2*sigma_const, color='g', linestyle='--')
plt.axhline(const+2*sigma_const, color='g', linestyle='--')
plt.grid()
plt.legend(loc='best', bbox_to_anchor=(0.5, -0.15, 0.5, 0.5))
plt.xlabel("Beam Energy [GeV]", fontsize=14)
plt.ylabel("Mean Radiation length for 10 electrons [cm]", fontsize=14)
plt.title("Mean Radiation length for electromagnetic showers\ncause by an electrons for various beam energy\nfrom simulated data with 95% CL", fontsize=14)
plt.savefig("Mean_radiation_length_electrons_G4_v2.svg")
