## THIS CODE IS IN PYTHON LANGUAGE ##

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# UTILITIES FUNCTION #

def lin_fit(x,a):
    return a*x
    
# DATA ACCESS AND PLOT #

data=[]
ne=[]
mean_data=[]
mean_ne=[]
mean_std_data=[]
mean_std_ne=[]

for i in range(1,6):
    i_file = np.loadtxt("charge_vs_ne_10e_"+str(i)+"GeV3_sigma_summed_charge_rectified2.txt", skiprows=1, delimiter=";")
    mean_data_elem=0
    mean_ne_elem=0
    temp_lst_sigma_data=[]
    temp_lst_sigma_ne=[]
    for j in range(10):
        data.append(i_file[:,0][j]*(1.6e-4)/(6.891e-3)) # IMPORTANT CORRECTION ON CONSTANTS
        ne.append(i_file[:,1][j]*1.6e-4)
        mean_data_elem+=i_file[:,0][j]*(1.6e-4)/(6.891e-3)
        mean_ne_elem+=i_file[:,1][j]*1.6e-4
        temp_lst_sigma_data.append(i_file[:,0][j]*(1.6e-4)/(6.891e-3))
        temp_lst_sigma_ne.append(i_file[:,1][j]*1.6e-4)
    print(mean_data_elem)
    mean_data.append(mean_data_elem/10)
    mean_ne.append(mean_ne_elem/10)

data = np.array(data)
ne = np.array(ne)
mean_data = np.array(mean_data)
mean_ne = np.array(mean_ne)




popt, cormat = curve_fit(lin_fit, mean_ne, mean_data, 1)

plt.plot(ne, data, "k+", label="Data")
plt.plot(ne, lin_fit(ne, popt), "-r", label=f"Best linear fit with p={popt[0]:.2f}")
plt.xlabel("Total deposited charge [fC]", fontsize=14)
plt.ylabel("Total reconstructed charge [fC]", fontsize=14)
plt.grid()
plt.legend()
plt.title("Total deposited charge VS the total reconstructed charge", fontsize=14)
plt.savefig("qtot_vs_ne_summed_charge.svg")
plt.savefig("qtot_vs_ne_summed_charge.pdf")
