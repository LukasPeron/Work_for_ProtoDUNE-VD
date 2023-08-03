### THIS CODE IS IN PYTHON LANGUAGE ###

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

## PLOT DATAS FOR THETA = 40° AND VARIABLE PHI

d_fin = []
d_event_fin = []
theta_fin = []
phi_fin = []
sigma_fin = []
track_length_fin = []
track_number_fin = []
d_for_track_length_fin = []

angle = 0
for i in range(6):
    track_length_data = np.loadtxt("theta_90_phi_var/2muon_10GeV_track_lenght_"+str(angle)+".txt", skiprows=1, delimiter=";")
    track_length = track_length_data[:,1]
    track_length_event = track_length_data[:,0]
    for h in range(len(track_length)):
        track_length_fin.append(track_length[h])

    track_number_data = np.loadtxt("theta_90_phi_var/2muon_10GeV_track_size_"+str(angle)+".txt", skiprows=1, delimiter=";")
    track_number = track_number_data[:,1]
    track_number_event = track_number_data[:,0]
    for l in range(len(track_number)):
        track_number_fin.append(track_number[l])

    data = np.loadtxt("theta_90_phi_var/2muon_10GeV_efficiency_"+str(angle)+".txt", skiprows=1, delimiter=";")
    d = data[:,0]
    sigma = data[:,1]
    theta = data[:,2]
    phi = data[:,3]
    d_event = data[:,4]
    for j in range(len(d)):
        d_fin.append(d[j])
        theta_fin.append(theta[j])
        phi_fin.append(phi[j])
        sigma_fin.append(sigma[j])
        d_event_fin.append(d_event[j])

    for event in range(10):
        for loop in range(int(track_number[event])):
            d_for_track_length_fin.append(data[event,0])
    angle+=15

def const_fit(x,a):
    const_array = np.linspace(np.min(x), np.max(x), np.size(x))
    for i in range(np.size(x)):
        const_array[i]=a
    return const_array

popt_const1, cormat_const1 = curve_fit(const_fit, theta_fin,d_fin, 10)
popt_const2, cormat_const2 = curve_fit(const_fit, phi_fin,d_fin, 10)


const1 = const_fit(theta_fin, *popt_const1)[0]
const2 = const_fit(phi_fin, *popt_const2)[0]

sigma_const1 = np.sqrt(cormat_const1[0,0])
sigma_const2 = np.sqrt(cormat_const2[0,0])

cumul_track = []
for i in range(1, int(np.max(track_number_fin))+1):
    elem=0
    for j in track_number_fin:
        if int(j)==i:
            elem+=i
    cumul_track.append(float(elem))
print(len(cumul_track))
print(np.sum(cumul_track))
cumul_track=np.array(cumul_track)
cumul_track/=np.sum(track_number_fin)
x_test = np.linspace(1,11,11)

# plt.figure(1)
# ax = plt.axes(projection="3d")
# ax.plot3D(theta_fin, phi_fin, d_fin, "ko")
# ax.set_xlabel("$\\theta$ [degree]")
# ax.set_ylabel("$\phi$ [degree]")
# ax.set_zlabel("$<d>$ [cm]")
# plt.title("Mean distance of reconstructed pandora track from simulated track")

plt.figure(2)
plt.xlim(np.min(theta_fin)-10,np.max(theta_fin)+10)
plt.plot(theta_fin, d_fin, marker='o', color='k', linestyle='', label="Data")
plt.hlines(const1, np.min(theta_fin)-10, np.max(theta_fin)+10, 'red', '-', label=f"Best constant fit $<d> = {popt_const1[0]:.2f}$ cm")
plt.fill_between([np.min(theta_fin)-10,np.max(theta_fin)+10], const1-sigma_const1, const1+sigma_const1, facecolor='g', alpha=0.2)
plt.fill_between([np.min(theta_fin)-10,np.max(theta_fin)+10], const1-2*sigma_const1, const1+2*sigma_const1, facecolor='g', alpha=0.2)
plt.hlines(const1-2*sigma_const1, np.min(theta_fin)-10, np.max(theta_fin)+10, color='g', linestyle='--')
plt.hlines(const1+2*sigma_const1, np.min(theta_fin)-10, np.max(theta_fin)+10, color='g', linestyle='--')
plt.hlines(2.8, np.min(theta_fin)-10, np.max(theta_fin)+10, 'blue', 'dashdot', label="$x$ systematic shift = 2.8 cm")
plt.grid()
plt.xlabel("$\\theta$ [degree]", fontsize=14)
plt.ylabel("$\langle d \\rangle$ [cm]", fontsize=14)
plt.title("$\\theta$ versus $\langle d \\rangle$ for $\\theta_{xz}=50$° and variable $\\theta_{yz}$", fontsize=14)
plt.legend(loc="upper left")
plt.savefig("theta_90_phi_var/d_vs_theta_usual.pdf")
plt.savefig("theta_90_phi_var/d_vs_theta_usual.svg")

plt.figure(3)
plt.xlim(np.min(phi_fin)-10,np.max(phi_fin)+10)
plt.plot(phi_fin, d_fin, marker='o', color='k', linestyle='', label="Data")
plt.grid()
plt.xlabel("$\phi$ [degree]", fontsize=14)
plt.ylabel("$\langle d \\rangle$ [cm]", fontsize=14)
plt.title("$\phi$ versus $\langle d \\rangle$ for $\\theta_{xz}=50$° and variable $\\theta_{yz}$", fontsize=14)
plt.hlines(const2, np.min(phi_fin)-10, np.max(phi_fin)+10, 'red', '-', label=f"Best constant fit $<d> = {popt_const2[0]:.2f}$ cm")
plt.hlines(2.8, np.min(phi_fin)-10, np.max(phi_fin)+10, 'blue', 'dashdot', label="$x$ systematic shift = 2.8 cm")
plt.fill_between([np.min(phi_fin)-10,np.max(phi_fin)+10], const2-sigma_const2, const2+sigma_const2, facecolor='g', alpha=0.2)
plt.fill_between([np.min(phi_fin)-10,np.max(phi_fin)+10], const2-2*sigma_const2, const2+2*sigma_const2, facecolor='g', alpha=0.2)
plt.hlines(const2-2*sigma_const2, np.min(phi_fin)-10, np.max(phi_fin)+10, color='g', linestyle='--')
plt.hlines(const2+2*sigma_const2, np.min(phi_fin)-10, np.max(phi_fin)+10, color='g', linestyle='--')
plt.legend(loc="best")
plt.savefig("theta_90_phi_var/d_vs_phi_usual.pdf")
plt.savefig("theta_90_phi_var/d_vs_phi_usual.svg")

fig, ax1 = plt.subplots()
ax1.set_xlabel("Number of pandora tracks", fontsize=14)
ax1.set_ylabel("$\langle d \\rangle$ [cm]", fontsize=14)
ln1 = ax1.plot(track_number_fin, d_fin, marker="o", color="k", linestyle="", label="Data")
ax1.tick_params(axis="y")
ax1.hlines(const1, np.min(track_number_fin)-1, np.max(track_number_fin)+1, 'red', '-', label=f"Best constant fit $<d> = {popt_const1[0]:.2f}$ cm")
plt.fill_between([np.min(track_number_fin)-1,np.max(track_number_fin)+1], const1-sigma_const1, const1+sigma_const1, facecolor='g', alpha=0.2)
plt.fill_between([np.min(track_number_fin)-1,np.max(track_number_fin)+1], const1-2*sigma_const1, const1+2*sigma_const1, facecolor='g', alpha=0.2)
plt.hlines(const1-2*sigma_const1, np.min(track_number_fin)-1, np.max(track_number_fin)+1, color='g', linestyle='--')
plt.hlines(const1+2*sigma_const1, np.min(track_number_fin)-1, np.max(track_number_fin)+1, color='g', linestyle='--')
ax2 = ax1.twinx()
ax2.set_ylabel("Proportional number of tracks", fontsize=14, color="b")
ln3 = ax2.plot(x_test, cumul_track, "-.ob", label="Proportional number of tracks")
ax2.tick_params(axis="y")
ax1.set_ylim(1, 20)
ax1.set_xlim(np.min(track_number_fin)-1, np.max(track_number_fin)+1)
ax1.legend(loc="best")
ax1.grid()
plt.title("Number of pandora tracks versus $\langle d \\rangle$ for\n$\\theta_{xz}=50$° and variable $\\theta_{yz}$", fontsize=14)
plt.show()
fig.savefig("theta_90_phi_var/d_vs_NTrack_usual.svg")
fig.savefig("theta_90_phi_var/d_vs_NTrack_usual.pdf")

plt.figure(5)
plt.ylim(0,70)
plt.xlim(np.min(track_length_fin)-10, np.max(track_length_fin)+10)
plt.plot(track_length_fin, d_for_track_length_fin, "ko", label="Datas")
plt.hlines(const1, np.min(track_length_fin)-10, np.max(track_length_fin)+10, 'red', '-', label=f"Best constant fit $<d> = {popt_const1[0]:.2f}$ cm")
plt.fill_between([np.min(track_length_fin)-10,np.max(track_length_fin)+10], const1-sigma_const1, const1+sigma_const1, facecolor='g', alpha=0.2)
plt.fill_between([np.min(track_length_fin)-10,np.max(track_length_fin)+10], const1-2*sigma_const1, const1+2*sigma_const1, facecolor='g', alpha=0.2)
plt.hlines(const1-2*sigma_const1, np.min(track_length_fin)-10, np.max(track_length_fin)+10, color='g', linestyle='--')
plt.hlines(const1+2*sigma_const1, np.min(track_length_fin)-10, np.max(track_length_fin)+10, color='g', linestyle='--')
plt.xlabel("Track length [cm]", fontsize=14)
plt.ylabel("$\langle d \\rangle$ [cm]", fontsize=14)
plt.title("Track length versus $\langle d \\rangle$ for $\\theta=40$° and variable $\phi$", fontsize=14)
plt.legend()
plt.grid()
plt.savefig("theta_90_phi_var/d_vs_TrackLength_usual.pdf")
plt.savefig("theta_90_phi_var/d_vs_TrackLength_usual.svg")

## PLOT DATAS FOR THETA = 90° AND VARIABLE PHI

d_fin = []
d_event_fin = []
theta_fin = []
phi_fin = []
sigma_fin = []
track_length_fin = []
track_number_fin = []
d_for_track_length_fin = []

def isnan(num):
    return num!=num

angle = 0
for i in range(6):
    data = np.loadtxt("theta_90_phi_var/10muon_10GeV_efficiency_"+str(angle)+".txt", skiprows=1, delimiter=";")
    d = data[:,0]
    sigma = data[:,1]
    theta = data[:,2]
    phi = data[:,3]
    d_event = data[:,4]
    for j in range(len(d)):
        if isnan(d[j]):
            continue
        if theta[j]<80 or theta[j]>100:
            continue
        else:
            d_fin.append(d[j])
            theta_fin.append(theta[j])
            phi_fin.append(phi[j])
            sigma_fin.append(sigma[j])
            d_event_fin.append(d_event[j])

    track_length_data = np.loadtxt("theta_90_phi_var/10muon_10GeV_track_lenght_"+str(angle)+".txt", skiprows=1, delimiter=";")
    track_length = track_length_data[:,1]
    track_length_event = track_length_data[:,0]
    for h in range(len(track_length)):
        event = int(track_length_event[h])
        theta_value = theta[event]
        if theta_value<80 or theta_value>100:
            continue
        else:
            if isnan(track_length[h]):
                continue
            else:
                track_length_fin.append(track_length[h])

    track_number_data = np.loadtxt("theta_90_phi_var/10muon_10GeV_track_size_"+str(angle)+".txt", skiprows=1, delimiter=";")
    track_number = track_number_data[:,1]
    track_number_event = track_number_data[:,0]
    for l in range(len(track_number)):
        if track_number[l]==0:
            continue
        if theta[l]<80 or theta[l]>100:
            continue
        else:            
            track_number_fin.append(track_number[l])

    # track_number_data = np.loadtxt("theta_90_phi_var/10muon_10GeV_track_size_"+str(angle)+".txt", skiprows=1, delimiter=";")
    # track_number = track_number_data[:,1]
    # track_number_event = track_number_data[:,0]
    # for l in range(len(track_number)):
    #     if track_number[l]==0:
    #         continue
    #     else:            
    #         track_number_fin.append(track_number[l])

    for event in range(len(track_number)):
        if theta[event]<80 or theta[event]>100:
            continue
        else:
            for loop in range(int(track_number[event])):
                d_for_track_length_fin.append(data[event,0])
    
    angle+=15

def const_fit(x,a):
    const_array = np.linspace(np.min(x), np.max(x), np.size(x))
    for i in range(np.size(x)):
        const_array[i]=a
    return const_array

popt_const1, cormat_const1 = curve_fit(const_fit, theta_fin,d_fin, 10)
popt_const2, cormat_const2 = curve_fit(const_fit, phi_fin,d_fin, 10)


const1 = const_fit(theta_fin, *popt_const1)[0]
const2 = const_fit(phi_fin, *popt_const2)[0]

sigma_const1 = np.sqrt(cormat_const1[0,0])
sigma_const2 = np.sqrt(cormat_const2[0,0])

cumul_track = []
for i in range(1, int(np.max(track_number_fin))+1):
    elem=0
    for j in track_number_fin:
        if int(j)==i:
            elem+=i
    cumul_track.append(float(elem))
# print(len(cumul_track))
# print(np.sum(cumul_track))
cumul_track=np.array(cumul_track)
cumul_track/=np.sum(track_number_fin)
x_test = np.linspace(1,int(np.max(track_number_fin)),int(np.max(track_number_fin)))

# plt.figure(1)
# ax = plt.axes(projection="3d")
# ax.plot3D(theta_fin, phi_fin, d_fin, "ko")
# ax.set_xlabel("$\\theta$ [degree]")
# ax.set_ylabel("$\phi$ [degree]")
# ax.set_zlabel("$<d>$ [cm]")
# plt.title("Mean distance of reconstructed pandora track from simulated track")

plt.figure(2)
plt.xlim(np.min(theta_fin)-10,np.max(theta_fin)+10)
plt.plot(theta_fin, d_fin, marker='o', color='k', linestyle='', label="Data")
plt.hlines(const1, np.min(theta_fin)-10, np.max(theta_fin)+10, 'red', '-', label=f"Best constant fit $<d> = {popt_const1[0]:.2f}$ cm")
plt.fill_between([np.min(theta_fin)-10,np.max(theta_fin)+10], const1-sigma_const1, const1+sigma_const1, facecolor='g', alpha=0.2)
plt.fill_between([np.min(theta_fin)-10,np.max(theta_fin)+10], const1-2*sigma_const1, const1+2*sigma_const1, facecolor='g', alpha=0.2)
plt.hlines(const1-2*sigma_const1, np.min(theta_fin)-10, np.max(theta_fin)+10, color='g', linestyle='--')
plt.hlines(const1+2*sigma_const1, np.min(theta_fin)-10, np.max(theta_fin)+10, color='g', linestyle='--')
plt.hlines(2.8, np.min(theta_fin)-10, np.max(theta_fin)+10, 'blue', 'dashdot', label="$x$ systematic shift = 2.8 cm")
plt.grid()
plt.xlabel("$\\theta$ [degree]", fontsize=14)
plt.ylabel("$\langle d \\rangle$ [cm]", fontsize=14)
plt.title("$\\theta$ versus $\langle d \\rangle$ for $\\theta_{xz}=0$° and variable $\\theta_{yz}$", fontsize=14)
plt.legend(loc=l, bbox_to_anchor=(0.3,0.8))
plt.savefig("theta_90_phi_var/d_vs_theta_horizontal.pdf")
plt.savefig("theta_90_phi_var/d_vs_theta_horizontal.svg")

plt.figure(3)
plt.xlim(np.min(phi_fin)-10,np.max(phi_fin)+10)
plt.plot(phi_fin, d_fin, marker='o', color='k', linestyle='', label="Data")
plt.grid()
plt.xlabel("$\phi$ [degree]", fontsize=14)
plt.ylabel("$\langle d \\rangle$ [cm]", fontsize=14)
plt.title("$\phi$ versus $\langle d \\rangle$ for $\\theta_{xz}=0$° and variable $\\theta_{yz}$", fontsize=14)
plt.hlines(const2, np.min(phi_fin)-10, np.max(phi_fin)+10, 'red', '-', label=f"Best constant fit $<d> = {popt_const2[0]:.2f}$ cm")
plt.hlines(2.8, np.min(phi_fin)-10, np.max(phi_fin)+10, 'blue', 'dashdot', label="$x$ systematic shift = 2.8 cm")
plt.fill_between([np.min(phi_fin)-10,np.max(phi_fin)+10], const2-sigma_const2, const2+sigma_const2, facecolor='g', alpha=0.2)
plt.fill_between([np.min(phi_fin)-10,np.max(phi_fin)+10], const2-2*sigma_const2, const2+2*sigma_const2, facecolor='g', alpha=0.2)
plt.hlines(const2-2*sigma_const2, np.min(phi_fin)-10, np.max(phi_fin)+10, color='g', linestyle='--')
plt.hlines(const2+2*sigma_const2, np.min(phi_fin)-10, np.max(phi_fin)+10, color='g', linestyle='--')
plt.legend(loc="best")
plt.savefig("theta_90_phi_var/d_vs_phi_horizontal.pdf")
plt.savefig("theta_90_phi_var/d_vs_phi_horizontal.svg")

fig, ax1 = plt.subplots()
ax1.set_xlabel("Number of pandora tracks", fontsize=14)
ax1.set_ylabel("$\langle d \\rangle$ [cm]", fontsize=14)
ln1 = ax1.plot(track_number_fin, d_fin, marker="o", color="k", linestyle="", label="Data")
ax1.tick_params(axis="y")
ax1.hlines(const1, np.min(track_number_fin)-1, np.max(track_number_fin)+1, 'red', '-', label=f"Best constant fit $<d> = {popt_const1[0]:.2f}$ cm")
ax1.fill_between([np.min(track_number_fin)-1,np.max(track_number_fin)+1], const1-sigma_const1, const1+sigma_const1, facecolor='g', alpha=0.2)
ax1.fill_between([np.min(track_number_fin)-1,np.max(track_number_fin)+1], const1-2*sigma_const1, const1+2*sigma_const1, facecolor='g', alpha=0.2)
ax1.hlines(const1-2*sigma_const1, np.min(track_number_fin)-1, np.max(track_number_fin)+1, color='g', linestyle='--')
ax1.hlines(const1+2*sigma_const1, np.min(track_number_fin)-1, np.max(track_number_fin)+1, color='g', linestyle='--')
ax2 = ax1.twinx()
ax2.set_ylabel("Proportional number of tracks", fontsize=14, color="b")
ln3 = ax2.plot(x_test, cumul_track, "-.ob", label="Proportional number of tracks")
ax2.tick_params(axis="y")
ax1.set_ylim(1, 20)
ax1.set_xlim(np.min(track_number_fin)-1, np.max(track_number_fin)+1)
ax1.legend(loc="upper right")
ax1.grid()
plt.title("Number of pandora tracks versus $\langle d \\rangle$ for\n$\\theta_{xz}=0$° and variable $\\theta_{yz}$", fontsize=14)
plt.show()
fig.savefig("theta_90_phi_var/d_vs_NTrack_horizontal.svg")
fig.savefig("theta_90_phi_var/d_vs_NTrack_horizontal.pdf")

plt.figure(5)
plt.ylim(0,70)
plt.xlim(np.min(track_length_fin)-10, np.max(track_length_fin)+10)
plt.plot(track_length_fin, d_for_track_length_fin, "ko", label="Datas")
plt.hlines(const1, np.min(track_length_fin)-10, np.max(track_length_fin)+10, 'red', '-', label=f"Best constant fit $<d> = {popt_const1[0]:.2f}$ cm")
plt.fill_between([np.min(track_length_fin)-10,np.max(track_length_fin)+10], const1-sigma_const1, const1+sigma_const1, facecolor='g', alpha=0.2)
plt.fill_between([np.min(track_length_fin)-10,np.max(track_length_fin)+10], const1-2*sigma_const1, const1+2*sigma_const1, facecolor='g', alpha=0.2)
plt.hlines(const1-2*sigma_const1, np.min(track_length_fin)-10, np.max(track_length_fin)+10, color='g', linestyle='--')
plt.hlines(const1+2*sigma_const1, np.min(track_length_fin)-10, np.max(track_length_fin)+10, color='g', linestyle='--')
plt.xlabel("Track length [cm]", fontsize=14)
plt.ylabel("$\langle d \\rangle$ [cm]", fontsize=14)
plt.title("Track length versus $\langle d \\rangle$ for $\\theta=90$° and variable $\phi$", fontsize=14)
plt.legend()
plt.grid()
plt.savefig("theta_90_phi_var/d_vs_TrackLength_horizontal.pdf")
plt.savefig("theta_90_phi_var/d_vs_TrackLength_horizontal.svg")
