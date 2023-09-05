### THIS CODE IS IN PYTHON LANGUAGE ###

import matplotlib.pyplot as plt
import numpy as np

# plot de nb_hits et q_tot en fonction de sigma_treshold et E_electron

nb_hits_5GeV_intgl = np.array([1418, 1307])
q_tot_5GeV_intgl = np.array([500303, 494615])

nb_hits_1GeV_intgl = np.array([386, 362, 310])
q_tot_1GeV_intgl = np.array([91049.9, 89676.2, 84769.1])

nb_hits_5GeV_SADC = np.array([450, 418])
q_tot_5GeV_SADC = np.array([505652, 502068])

nb_hits_1GeV_SADC = np.array([194, 185, 171])
q_tot_1GeV_SADC = np.array([99610.6, 98978.8, 95826.4])


nb_sigma_5GeV = np.array([3, 5])
nb_sigma_1GeV = np.array([2, 3, 5])

fig2, ax3 = plt.subplots() 
# plt.title("Total charge in the reconstruction versus \n the reconstruction treshold for 1 GeV electron", fontsize=14)
ax3.set_xlabel('Reconstruction trehsold [$\sigma$]', fontsize=14) 
ax3.set_ylabel('Total charge in the reconstruction [ADC]', fontsize=14) 
ax3.plot(nb_sigma_1GeV, q_tot_1GeV_intgl, 'ro--', label="Integral method")
ax3.plot(nb_sigma_1GeV, q_tot_1GeV_SADC, 'bo--', label="SummedADC method")
plt.grid()
plt.legend()
plt.show()
fig2.savefig("1GeV_q_tot.svg")
fig2.savefig("1GeV_q_tot.pdf")
