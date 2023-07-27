## THIS CODE IS IN PYTHON LANGUAGE ##

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

data=[]
ne=[]

for i in range(1,6):
    i_file = np.loadtxt("charge_vs_ne_10e_"+str(i)+"GeV3_sigma_summed_charge_rectified2.txt", skiprows=1, delimiter=";")
    for j in range(10):
        data.append(i_file[:,0][j]*(1.6e-4)/(6.891e-3)) # CORRECTION IMPORTANTE SUR LES CONSTANTES
        ne.append(i_file[:,1][j]*1.6e-4)

data = np.sort(np.array(data))
ne = np.sort(np.array(ne))


def lin_fit(x,a):
    return a*x


def aff_fit(x,a,b):
    return a*x + b

popt, cormat = curve_fit(lin_fit, ne, data, 1)
# popt2, cormat2 = curve_fit(aff_fit, ne, data, (1,0))

plt.plot(ne, data, "k+", label="data")
plt.plot(ne, lin_fit(ne, popt), "-r", label=f"best linear fit with p={popt[0]:.2f}")
# plt.plot(ne, aff_fit(ne, *popt2), "-b", label="best affine fit")
plt.xlabel("Total deposited charge [fC]", fontsize=14)
plt.ylabel("Total reconstructed charge [fC]", fontsize=14)
plt.grid()
plt.legend()
plt.title("Total deposited charge VS the total reconstructed charge", fontsize=14)
plt.savefig("qtot_vs_ne_summed_charge.svg")
