### THIS CODE IS IN PYTHON LANGUAGE ###

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import numpy as np
from scipy.optimize import curve_fit

# UTILITIES CLASS AND FUNCTIONS #

class nlcmap(object):
    def __init__(self, cmap, levels):
        self.cmap = cmap
        self.N = cmap.N
        self.monochrome = self.cmap.monochrome
        self.levels = np.asarray(levels, dtype='float64')
        self._x = self.levels
        self.levmax = self.levels.max()
        self.transformed_levels = np.linspace(0.0, self.levmax,
             len(self.levels))

    def __call__(self, xi, alpha=1.0, **kw):
        yi = np.interp(xi, self._x, self.transformed_levels)
        return self.cmap(yi / self.levmax, alpha)

def const_fit(x,a):
    const_array = np.linspace(np.min(x), np.max(x), np.size(x))
    for i in range(np.size(x)):
        const_array[i]=a
    return const_array

def isnan(num):
    return num!=num

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
plt.title("$\\theta$ versus $\langle d \\rangle$ for Dataset 1", fontsize=14)
plt.legend(loc="upper left")
plt.savefig("theta_90_phi_var/d_vs_theta_usual.pdf")
plt.savefig("theta_90_phi_var/d_vs_theta_usual.svg")

plt.figure(3)
plt.xlim(np.min(phi_fin)-10,np.max(phi_fin)+10)
plt.plot(phi_fin, d_fin, marker='o', color='k', linestyle='', label="Data")
plt.grid()
plt.xlabel("$\phi$ [degree]", fontsize=14)
plt.ylabel("$\langle d \\rangle$ [cm]", fontsize=14)
plt.title("$\phi$ versus $\langle d \\rangle$ for Dataset 1", fontsize=14)
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
plt.title("Number of pandora tracks versus $\langle d \\rangle$ for Dataset 1", fontsize=14)
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
plt.title("Track length versus $\langle d \\rangle$ for Dataset 1", fontsize=14)
plt.legend()
plt.grid()
plt.savefig("theta_90_phi_var/d_vs_TrackLength_usual.pdf")
plt.savefig("theta_90_phi_var/d_vs_TrackLength_usual.svg")

track_length_per_event = []
seuil = 0
for nb_track in track_number_fin:
    elem_lst = []
    for i in range(seuil, seuil+int(nb_track)):
        elem_lst.append(track_length_fin[i])
        seuil+=1
    track_length_per_event.append(elem_lst)

plt.figure(6)
for elem in range(len(theta_fin)):
    # plt.plot(theta_fin[elem], [track_length_per_event[elem]], 'ok')
    plt.plot(theta_fin[elem], sum(x for x in track_length_per_event[elem]), "ok")
plt.grid()
# plt.legend(["Data", "Total track length"])
plt.xlabel("$\\theta$ [degree]", fontsize=14)
plt.ylabel("Track length [cm]", fontsize=14)

plt.figure(7)
for elem in range(len(theta_fin)):
    plt.plot(phi_fin[elem], [track_length_per_event[elem]], 'ok')
plt.grid()
plt.xlabel("$\phi$ [degree]", fontsize=14)
plt.ylabel("Track length [cm]", fontsize=14)

plt.figure(8)
plt.plot(theta_fin, track_number_fin, 'ko')
plt.grid()
plt.xlabel("$\\theta$ [degree]", fontsize=14)
plt.ylabel("Number of track", fontsize=14)

## PLOT DATAS FOR THETA = 90° AND VARIABLE PHI

d_fin = []
d_event_fin = []
theta_fin = []
phi_fin = []
sigma_fin = []
track_length_fin = []
track_number_fin_d2 = []
d_for_track_length_fin = []

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
        else:            
            track_number_fin_d2.append(track_number[l])

    for event in range(len(track_number)):
        for loop in range(int(track_number[event])):
            d_for_track_length_fin.append(data[event,0])
    
    angle+=15

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
x_test = np.linspace(1,int(np.max(track_number_fin)),int(np.max(track_number_fin)))


plt.figure(9)
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
plt.title("$\\theta$ versus $\langle d \\rangle$ for Dataset 2", fontsize=14)
plt.legend(loc=l, bbox_to_anchor=(0.3,0.8))
plt.savefig("theta_90_phi_var/d_vs_theta_horizontal.pdf")
plt.savefig("theta_90_phi_var/d_vs_theta_horizontal.svg")

plt.figure(10)
plt.xlim(np.min(phi_fin)-10,np.max(phi_fin)+10)
plt.plot(phi_fin, d_fin, marker='o', color='k', linestyle='', label="Data")
plt.grid()
plt.xlabel("$\phi$ [degree]", fontsize=14)
plt.ylabel("$\langle d \\rangle$ [cm]", fontsize=14)
plt.title("$\phi$ versus $\langle d \\rangle$ for Dataset 2", fontsize=14)
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
plt.title("Number of pandora tracks versus $\langle d \\rangle$ for Dataset 2", fontsize=14)
plt.show()
fig.savefig("theta_90_phi_var/d_vs_NTrack_horizontal.svg")
fig.savefig("theta_90_phi_var/d_vs_NTrack_horizontal.pdf")

plt.figure(11)
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
plt.title("Track length versus $\langle d \\rangle$ for Dataset 2", fontsize=14)
plt.legend()
plt.grid()
plt.savefig("theta_90_phi_var/d_vs_TrackLength_horizontal.pdf")
plt.savefig("theta_90_phi_var/d_vs_TrackLength_horizontal.svg")

track_length_per_event_d2 = []
seuil = 0
for nb_track in track_number_fin_d2:
    elem_lst = []
    for i in range(seuil, seuil+int(nb_track)):
        elem_lst.append(track_length_fin[i])
        seuil+=1
    track_length_per_event_d2.append(elem_lst)

plt.figure(12)
for elem in range(len(theta_fin)):
    # plt.plot(theta_fin[elem], [track_length_per_event[elem]], 'ok')
    plt.plot(theta_fin[elem], sum(x for x in track_length_per_event[elem]), "ko")
plt.grid()
# plt.legend(["Data", "Total track length"])
plt.xlabel("$\\theta$ [degree]", fontsize=14)
plt.ylabel("Track length [cm]", fontsize=14)

plt.figure(13)
for elem in range(len(theta_fin)):
    plt.plot(phi_fin[elem], [track_length_per_event[elem]], 'ok')
plt.grid()
plt.xlabel("$\phi$ [degree]", fontsize=14)
plt.ylabel("Track length [cm]", fontsize=14)

plt.figure(14)
plt.plot(theta_fin, track_number_fin, 'ko')
plt.grid()
plt.xlabel("$\\theta$ [degree]", fontsize=14)
plt.ylabel("Number of track", fontsize=14)

## TRACK LENGTH VS TRACK NUMBER ##

sorted_lst = [[], [], [], [], [], [], [], [], [], [], []]
sorted_lst_2 = [[], [], [], [], [], [], [], [], [], [], [], [], [], []]
mean_lst = []
mean_lst_2 = []
sigma_lst = []
sigma_lst_2 = []

for i in range(len(track_number_fin)):
    index = int(track_number_fin[i])-1
    elem_lst = np.mean(track_length_per_event[i])
    sorted_lst[index].append(elem_lst)

for i in range(11):
    elem = sorted_lst[i]
    mean = np.mean(elem)
    sigma = np.std(elem)
    mean_lst.append(mean)
    sigma_lst.append(sigma)

for i in range(len(track_number_fin_d2)):
    index = int(track_number_fin_d2[i])-1
    elem_lst = np.mean(track_length_per_event_d2[i])
    sorted_lst_2[index].append(elem_lst)

for i in range(14):
    elem = sorted_lst_2[i]
    if elem==[]:
        continue
    else:
        mean = np.mean(elem)
        sigma = np.std(elem)
    mean_lst_2.append(mean)
    sigma_lst_2.append(sigma)

x = np.linspace(1,11,11,endpoint=True)
x2 = np.linspace(1,14,14,endpoint=True)

prod = x*np.array(mean_lst)
prod_2 = np.array(np.linspace(1,len(mean_lst_2),len(mean_lst_2),endpoint=True)*np.array(mean_lst_2))

popt1, cormat1 = curve_fit(const_fit, x, prod, 400)
popt2, cormat2 = curve_fit(const_fit, np.linspace(1,8,8,endpoint=True), prod_2, 400)

const_2 = const_fit(np.linspace(1,8,8,endpoint=True), *popt2)
const_1 = const_fit(x, *popt1)

prod_2 = np.insert(prod_2, 6, None)
const_2 = np.insert(const_2, 6, None)
mean_lst_2 = np.insert(mean_lst_2, 6, None)
sigma_lst_2 = np.insert(sigma_lst_2, 6, None)
for i in range(5):
    prod_2 = np.insert(prod_2, 8, None)
    const_2 = np.insert(const_2, 8, None)
    mean_lst_2 = np.insert(mean_lst_2, 8, None)
    sigma_lst_2 = np.insert(sigma_lst_2, 8, None)

fig, ax = plt.subplots(1,1)


ax.errorbar(x, mean_lst, yerr=sigma_lst, marker='o', color='b', linestyle="", label="Dataset 1")
ax.errorbar(x2, mean_lst_2, yerr=sigma_lst_2, marker='o', color='r', linestyle="", alpha=0.8, label="Dataset 2")
ax.plot(x, prod, '+b', label="$N_{tracks} \\times L_{tracks}$ Dataset 1")
ax.plot(x2, prod_2, '+r', label="$N_{tracks} \\times L_{tracks}$ Dataset 2")
ax.hlines(const_1[0], 0, 15, color="b", linestyle="--")#, label="Best constant fit $\sum L_{tracks}$ ="+f"{const_1[0]:.0f}cm")
ax.hlines(const_2[0], 0, 15, color="r", linestyle="--")#, label="Best constant fit $\sum L_{tracks}$ ="+ f"{const_2[0]:.0f}cm")
ax.text(10, 345, "Best constant fit\n$N_{tracks} \\times L_{tracks}$ ="+f"{const_1[0]:.0f}cm\nDataset 1", fontsize=10, color="b")
ax.text(10, 205, "Best constant fit\n$N_{tracks} \\times L_{tracks}$ ="+f"{const_2[0]:.0f}cm\nDataset 2", fontsize=10, color="r")
ax.grid()
ax.set_xlabel("Number of pandora tracks ($N_{tracks}$)", fontsize=14)
ax.set_ylabel("Mean tracks length ($L_{tracks}$) [cm]", fontsize=14)
ax.legend(loc=l, bbox_to_anchor=(0.65,0.44))
fig.savefig("track_length_vs_number.svg")

## <d> vs theta and phi ##
d_fin = []
d_event_fin = []
theta_fin = []
phi_fin = []
sigma_fin = []

angle = 0
for i in range(6):
    data = np.loadtxt("theta_90_phi_var/2muon_10GeV_efficiency_"+str(angle)+".txt", skiprows=1, delimiter=";")
    d = data[:,0]
    sigma = data[:,1]
    theta = data[:,2]
    phi = data[:,3]
    d_event = data[:,4]
    for j in range(len(d)):
        if isnan(d[j]):
            continue
        else:
            d_fin.append(d[j])
            theta_fin.append(theta[j])
            phi_fin.append(phi[j])
            sigma_fin.append(sigma[j])
            d_event_fin.append(d_event[j])
    angle+=15

print(len(d_fin))

angle=0
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
        else:
            d_fin.append(d[j])
            theta_fin.append(theta[j])
            phi_fin.append(phi[j])
            sigma_fin.append(sigma[j])
            d_event_fin.append(d_event[j])
    angle+=15

x = theta_fin
y = phi_fin
t1 = np.array(d_fin)
mean1 = np.mean(d_fin)
sigma1 = np.std(d_fin)
xlim = 80
ylim = 100
tmax = np.max(t1.max())
tmin = np.min(t1.min())
levels = np.concatenate(([20, tmax], np.linspace(tmin, mean1 + 1.75*sigma1, 6)))
levels = levels[levels <= tmax]
levels.sort()
N = 256
vals = np.ones((N, 4))
vals[:, 0] = np.linspace(0/256, 1, N)
vals[:, 1] = np.linspace(0/256, 0, N)
vals[:, 2] = np.linspace(256/256, 0, N)
vals[:, 3] = np.linspace(0.4, 1, N)
newcmp = ListedColormap(vals)

cmap_nonlin = nlcmap(plt.cm.YlOrRd, levels)
fig, ax1 = plt.subplots(1,1, figsize=(5.4,5))
cbar_ax = fig.add_axes([0.95, 0.1, 0.05, 0.78]) # [0.2, x>0, 0.6, 0.05] for save
sm = plt.cm.ScalarMappable(cmap=plt.cm.YlOrRd, norm=plt.Normalize(vmin=0, vmax=tmax))
sm._A = []
cbar = fig.colorbar(sm, cax=cbar_ax, location="right")
cbar.set_label(label="$\langle d \\rangle$ [cm]", fontsize=14)
cbar.set_ticks(cmap_nonlin.transformed_levels)
cbar.set_ticklabels(["%.2f" % lev for lev in levels])
ax1.scatter(x,y,edgecolors=cmap_nonlin(t1), s=15, linewidths=4)
ax1.set_xlim(np.min(theta_fin)-5, np.max(theta_fin)+30)
ax1.set_ylim(np.min(phi_fin)-5, np.max(phi_fin)+5)
ax1.grid()
ax1.set_ylabel("$\phi$ [degree]", fontsize=14)
ax1.set_xlabel("$\\theta$ [degree]", fontsize=14)
ax1.set_title("$\langle d \\rangle$ vs $\\theta$ and $\phi$", fontsize=14)
zax = ax1.inset_axes([0.635, 0.02, 0.35, 0.3])
zax.scatter(x,y,edgecolors=cmap_nonlin(t1), s=15, linewidths=4)
zax.set_xlim(80, 100)
zax.set_ylim(0, 100)
zax.grid()
xticks=zax.xaxis.get_major_ticks()
yticks=zax.yaxis.get_major_ticks()
for i in range(len(xticks)):
    xticks[i].label1.set_visible(False)
for i in range(len(yticks)):
    yticks[i].label1.set_visible(False)
ax1.indicate_inset_zoom(zax, edgecolor="k")
plt.savefig("theta_phi_d_rapport.svg")
