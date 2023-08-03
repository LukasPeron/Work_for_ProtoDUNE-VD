# CE CODE EST EN LANGAGE PYTHON #

import numpy as np
import matplotlib.pyplot as plt

a = np.array([[0.0830899, -0.08717743, 0.02610534, -2.74655e-3, 4.39504e-5, 9.05605e-6, -3.97621e-7],
                [0.265283, -0.10167009, 0.00701793, 2.371288e-3, -5.020251e-4, 3.6531e-5, -9.46044e-7],
                [2.18838e-3, -2.914205e-3, 1.26639e-3, -7.6598e-5, -1.58882e-5, 2.18716e-6, -7.49728e-8],
                [-4.48746e-5, 4.75329e-5, -1.43471e-5, 1.19661e-6, 5.7891e-8, -1.2617e-8, 4.633e-10],
                [6.29882e-7, -6.72311e-7, 2.62963e-7, -5.1862e-8, 5.692e-9, -3.29e-10, 7.7e-12]])

E_e = 1e9*np.array([1, 2, 3, 4, 5]) # eV
Z = 18
m = 6.62e-26 #kg
p = 1.44e3 #kg/m^3
V = 6**3 #m^3

sigma_lst=[]
E_gamma_lst = []
mfp_discret = []
md_discret = []

for l in range(5):
    E_gamma = E_e[l]*(1-1/np.e)
    E_gamma_lst.append(E_gamma)
    k=E_gamma/511e3
    sigma = 0
    for i in range(7):
        inter = 0
        for j in range(5):
            inter += a[j][i]*(Z**j)
        inter *= np.log(k)**i
        sigma += inter
    sigma_lst.append(sigma)
    mfp = 100*m/(p*sigma*1e-28)
    mfp_discret.append(mfp)
    md_discret.append(0.44*mfp + 0.56*100*0.14)

E_e_cont = np.linspace(1e9, 5e9, 10000)
E_gamma_cont_lst = []
sigma_lst_cont = []

for l in range(len(E_e_cont)):
    E_gamma_cont = E_e_cont[l]*(1-1/np.e)
    E_gamma_cont_lst.append(E_gamma_cont)
    k=E_gamma_cont/511e3
    sigma_cont = 0
    for i in range(7):
        inter = 0
        for j in range(5):
            inter += a[j][i]*(Z**j)
        inter *= np.log(k)**i
        sigma_cont += inter
    sigma_lst_cont.append(sigma_cont)

E_gamma_cont_lst=np.array(E_gamma_cont_lst)/1e9
E_gamma_lst = np.array(E_gamma_lst)/1e9

mean_free_path = np.zeros(len(sigma_lst_cont))
for i in range(len(mean_free_path)):
    mean_free_path[i] = 100*m/(p*sigma_lst_cont[i]*1e-28)
mean_distance = np.zeros(len(mean_free_path))
for i in range(len(mean_distance)):
    mean_distance[i] = 0.44*mean_free_path[i] + 0.56*100*0.14

plt.figure(1)
plt.plot(E_gamma_cont_lst/(1-1/np.e), sigma_lst_cont, "-r")
plt.plot(E_gamma_lst/(1-1/np.e), sigma_lst, "ok")
plt.xlabel("$E_{e^-}$ [GeV]", fontsize=14)
plt.ylabel("$\sigma$ [b/atom]", fontsize=14)
plt.title("$\sigma$ en fonction de l'énergie de l'électron primaire", fontsize=14)
plt.grid()
plt.savefig("sigma_vs_E_elec.svg")

plt.figure(2)
plt.plot(E_gamma_cont_lst/(1-1/np.e), mean_distance, '-r')
plt.plot(E_gamma_lst/(1-1/np.e), md_discret, "ok")
plt.xlabel("$E_{e^-}$ [GeV]", fontsize=14)
plt.ylabel("$\\bar{d_\gamma}$ [cm]", fontsize=14)
plt.title("$\\bar{d_\gamma}$ en fonction de l'énergie de l'électron primaire", fontsize=14)
plt.grid()
plt.savefig("Mean_distance_vs_E_elec.svg")
