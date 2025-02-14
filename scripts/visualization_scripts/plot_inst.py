# Python vizualization scripts, using Python 3.4.1
# Author: Jérémy Auffinger jeremy.auffinger.blackhawk@protonmail.com
# Last modification: 15 September 2023

## 1 - Necessary importations

import matplotlib.pyplot as plt
import numpy as np
#!python
import pylab
from pylab import arange,pi,sin,cos,sqrt
fig_width_pt = 170*2.83465  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height_pt = 160*2.83465
fig_height = fig_height_pt*inches_per_pt     # height in inches
fig_size =  [fig_width,fig_height]
params = {'backend': 'ps',
          'axes.labelsize': 10,
          'axes.linewidth': 0.5,
          'font.size': 10,
          'figure.titlesize': 10,
          'legend.fontsize': 9,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'text.usetex': True,
          'figure.figsize': fig_size}
pylab.rcParams.update(params)

## 2 - Folder definition

# Here put the BlackHawk path
BH_path = "[BlackHawk path]"
#
result_folder = BH_path + "/results/"
# Here put the name of the 'destination_folder'
data_folder = "[destination_folder]"
folder = result_folder + data_folder + "/"
# Here put a path to some figures folder
fig_folder = folder

## 3 - Recovering data

data_spectrum = np.genfromtxt(folder+"BH_spectrum.txt",skip_header = 2)
data_primary = np.genfromtxt(folder+"instantaneous_primary_spectra.txt",skip_header = 2)
data_secondary = np.genfromtxt(folder+"instantaneous_secondary_spectra.txt",skip_header = 2)

## 4 - Plotting BH spectrum

pylab.figure(1)
pylab.clf()
pylab.axes([0.125,0.2,0.90 - 0.125,0.90-0.2])
dens_min = np.min(data_spectrum[:,2])/2
dens_max = 2*np.max(data_spectrum[:,2])
widths_up = np.zeros(len(data_spectrum))
widths_down = np.zeros(len(data_spectrum))
for i in range(1,len(data_spectrum)-1):
    widths_down[i] = 0.5*(data_spectrum[i,0] - data_spectrum[i-1,0])
    widths_up[i] = 0.5*(data_spectrum[i+1,0] - data_spectrum[i,0])
widths_up[0] = 0.5*(data_spectrum[1,0] - data_spectrum[0,0])
widths_down[0] = 0.5*(data_spectrum[0,0])
widths_down[-1] = 0.5*(data_spectrum[-1,0] - data_spectrum[-2,0])
widths_up[-1] = 0.5*(data_spectrum[-1,0] - data_spectrum[-2,0])
widths = np.array([widths_down,widths_up])
plt.ylim(dens_min,dens_max)
plt.xlim(data_spectrum[0,0]/10,data_spectrum[-1,0]*10)
pylab.scatter(data_spectrum[:,0],data_spectrum[:,2],s=1,c='black')
pylab.yscale('log')
pylab.xscale('log')
pylab.errorbar(data_spectrum[:,0],data_spectrum[:,2],xerr = widths,yerr = None,elinewidth = 2,linewidth = 0)
pylab.xlabel('$M{\\rm \,\,(g)}$')
pylab.ylabel('${\\rm d}n\,\, ({\\rm cm}^{-3})$')
title = '${\\rm BH\,\, mass\,\,distribution\,\,}-{\\rm \,\,' + data_folder + '}$'
pylab.title(title,fontsize = 10,y = 1.02)
pylab.grid()

fig_name = fig_folder + data_folder + "_" + "spectrum.pdf"

pylab.savefig(fig_name)

## 5 - Plotting options (primary data)

# Put 1 to plot the particle spectrum
photon_primary=1
gluons_primary=0
higgs_primary=0
W_primary=0
Z_primary=0
neutrinos_primary=0
electron_primary=1
muon_primary=0
tau_primary=0
up_primary=0
down_primary=0
charm_primary=0
strange_primary=0
top_primary=0
bottom_primary=0
graviton_primary=0

part_show_primary=np.zeros(16)
part_show_primary=np.array(part_show_primary,int)
part_show_primary[0] = photon_primary
part_show_primary[1] = gluons_primary
part_show_primary[2] = higgs_primary
part_show_primary[3] = W_primary
part_show_primary[4] = Z_primary
part_show_primary[5] = neutrinos_primary
part_show_primary[6] = electron_primary
part_show_primary[7] = muon_primary
part_show_primary[8] = tau_primary
part_show_primary[9] = up_primary
part_show_primary[10] = down_primary
part_show_primary[11] = charm_primary
part_show_primary[12] = strange_primary
part_show_primary[13] = top_primary
part_show_primary[14] = bottom_primary
part_show_primary[15] = graviton_primary

labels_primary=np.array(np.zeros(16),str)
labels_primary[0]="$\gamma$"
labels_primary[1]="$g$"
labels_primary[2]="$h$"
labels_primary[3]="$W^\pm$"
labels_primary[4]="$Z^0$"
labels_primary[5]="$\\nu,\overline{\\nu}$"
labels_primary[6]="$e^\pm$"
labels_primary[7]="$\mu^\pm$"
labels_primary[8]="$\\tau^\pm$"
labels_primary[9]="$u,\overline{u}$"
labels_primary[10]="$d,\overline{d}$"
labels_primary[11]="$c,\overline{c}$"
labels_primary[12]="s,\overline{s}"
labels_primary[13]="t,\overline{t}"
labels_primary[14]="b,\overline{b}"
labels_primary[15]="${\\rm G}$"

## 6 - Plotting primary spectra

pylab.figure(2)
pylab.clf()
pylab.axes([0.125,0.2,0.90 - 0.125,0.90-0.2])
flux_max = 0.
for i in range(16):
    if part_show_primary[i]:
        flux_max = max(flux_max,max(data_primary[:,i+1]))
plt.ylim(flux_max/1e+10,flux_max*10.)
for i in range(16):
    if part_show_primary[i]:
        pylab.loglog(data_primary[:,0],data_primary[:,i+1],label = labels_primary[i],linewidth = 2)
pylab.xlabel('$E{\\rm \,\, (GeV)}$')
pylab.ylabel('${\\rm d}^2n/{\\rm d}t{\\rm d}E\,\, ({\\rm GeV}^{-1}\cdot{\\rm s}^{-1}\cdot{\\rm cm}^{-3})$')
title = '${\\rm Primary\,\, spectra\,\,}-{\\rm \,\,' + data_folder + '}$'
pylab.title(title,fontsize = 10,y = 1.02)
pylab.grid()
pylab.legend(loc = 1)

fig_name = fig_folder + data_folder + "_" + "primary.pdf"

pylab.savefig(fig_name)

## 7 - Plotting options (secondary data)

epoch = 1 # 0: BBN epoch, 1: present epoch

if epoch == 0:
    nb_fin_part = 11

    # Put 1 to plot the particle spectrum
    photon_secondary=1
    electron_secondary=1
    muon_secondary=0
    nu_e_secondary=1
    nu_mu_secondary=0
    nu_tau_secondary=0
    pipm_secondary=0
    K0L_secondary=0
    Kpm_secondary=0
    proton_secondary=0
    neutron_secondary=0

    part_show_secondary=np.zeros(nb_fin_part)
    part_show_secondary=np.array(part_show_secondary,int)
    part_show_secondary[0] = photon_secondary
    part_show_secondary[1] = electron_secondary
    part_show_secondary[2] = muon_secondary
    part_show_secondary[3] = nu_e_secondary
    part_show_secondary[4] = nu_mu_secondary
    part_show_secondary[5] = nu_tau_secondary
    part_show_secondary[6] = pipm_secondary
    part_show_secondary[7] = K0L_secondary
    part_show_secondary[8] = Kpm_secondary
    part_show_secondary[9] = proton_secondary
    part_show_secondary[10] = neutron_secondary

    labels_secondary=np.array(np.zeros(nb_fin_part),str)
    labels_secondary[0]="$\gamma$"
    labels_secondary[1]="$e^\pm$"
    labels_secondary[2]="$\mu^\pm$"
    labels_secondary[3]="$\\nu_e,\overline{\\nu}_e$"
    labels_secondary[4]="$\\nu_\mu,\overline{\\nu}_\mu$"
    labels_secondary[5]="$\\nu_\\tau,\overline{\\nu}_\\tau$"
    labels_secondary[6]="$\pi^\pm$"
    labels_secondary[7]="$K_{\\rm L}^0$"
    labels_secondary[8]="$K^\pm$"
    labels_secondary[9]="${\\rm p},\overline{\\rm p}$"
    labels_secondary[10]="${\\rm n},\overline{\\rm n}$"

elif epoch == 1:
    nb_fin_part = 6

    # Put 1 to plot the particle spectrum
    photon_secondary=1
    electron_secondary=1
    nu_e_secondary=0
    nu_mu_secondary=0
    nu_tau_secondary=0
    proton_secondary=1

    part_show_secondary=np.zeros(nb_fin_part)
    part_show_secondary=np.array(part_show_secondary,int)
    part_show_secondary[0] = photon_secondary
    part_show_secondary[1] = electron_secondary
    part_show_secondary[2] = nu_e_secondary
    part_show_secondary[3] = nu_mu_secondary
    part_show_secondary[4] = nu_tau_secondary
    part_show_secondary[5] = proton_secondary

    labels_secondary=np.array(np.zeros(nb_fin_part),str)
    labels_secondary[0]="$\gamma$"
    labels_secondary[1]="$e^\pm$"
    labels_secondary[2]="$\\nu_e,\overline{\\nu}_e$"
    labels_secondary[3]="$\\nu_\mu,\overline{\\nu}_\mu$"
    labels_secondary[4]="$\\nu_\\tau,\overline{\\nu}_\\tau$"
    labels_secondary[5]="${\\rm p},\overline{\\rm p}$"

## 8 - Plotting secondary spectra

pylab.figure(3)
pylab.clf()
pylab.axes([0.125,0.2,0.90 - 0.125,0.90-0.2])
flux_max = 0.
for i in range(nb_fin_part):
    if part_show_secondary[i]:
        flux_max = max(flux_max,max(data_secondary[:,i+1]))
plt.ylim(flux_max/1e+10,flux_max*10.)
for i in range(nb_fin_part):
    if part_show_secondary[i]:
        pylab.loglog(data_secondary[:,0],data_secondary[:,i+1],label = labels_secondary[i],linewidth = 2)
pylab.xlabel('$E{\\rm \,\, (GeV)}$')
pylab.ylabel('${\\rm d}^2n/{\\rm d}t{\\rm d}E\,\, ({\\rm GeV}^{-1}\cdot{\\rm s}^{-1}\cdot{\\rm cm}^{-3})$')
title = '${\\rm Secondary\,\, spectra\,\,}-{\\rm \,\,' + data_folder + '}$'
pylab.title(title,fontsize = 10,y = 1.02)
pylab.grid()
pylab.legend(loc = 1)

fig_name = fig_folder + data_folder + "_" + "secondary.pdf"

pylab.savefig(fig_name)

## 9 - Closing all figures
pylab.close("all")










