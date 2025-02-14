# Script written by J. Auffinger jeremy.auffinger.blackhawk@protonmail.com & A. Arbey alexandre.arbey@ens-lyon.fr
# Last modification 07 April 2024
## Hazma importations

import numpy as np
import hazma.decay
import hazma.decay_helper_functions
from hazma.decay import neutral_pion as dnde_pi0
from hazma.decay import muon as dnde_mu
from hazma.decay import charged_pion as dnde_pi
from hazma.positron_spectra import muon as dnde_p_mu
from hazma.positron_spectra import charged_pion as dnde_p_pi

for i in range(len(e_init)):
    table_mu_e[i,:] = dnde_p_mu(e_fin,e_init[i])
    table_pipm_e[i,:] = dnde_p_pi(e_fin,e_init[i])

file = "hamza/Hazma-master/Hazma-master/electron.txt"
print("%15s%15s%15s%15s%15s\n"%("initial_energy","final_energy","pi0>e+/-","pi+/->e+/-","mu+/->e+/-"),file=open(file,"w"),end = "")
for i in range(Enum_init):
    for j in range(Enum_fin):
        print("%15.5e%15.5e%15.5e%15.5e%15.5e\n"%(e_init[i]/1000.,e_fin[j]/1000.,table_pi0_e[i,j]*1000.,table_pipm_e[i,j]*1000.,table_mu_e[i,j]*1000.),file=open(file,"a"),end = "")
## Defining the quantities

Enum_init = 250 # number of initial energies
Enum_fin = 500 # number of final energies

Emin_init = 1e-6 # minimal initial energy in GeV
Emax_init = 5.0 # maximal initial energy in GeV
Emin_fin = 1e-6 # minimal final energy in GeV
Emax_fin = 5. # maximal final energy in GeV

e_init = np.logspace(np.log10(Emin_init),np.log10(Emax_init),Enum_init)*1000. # Hazma works in MeV
e_fin = np.logspace(np.log10(Emin_fin),np.log10(Emax_fin),Enum_fin)*1000. # Hazma works in MeV

table_pi0 = np.zeros([Enum_init,Enum_fin]) # contains the pi0 -> gamma decay rates dN/dE in MeV^(-1)
table_pipm = np.zeros([Enum_init,Enum_fin]) # contains the pi+- -> gamma decay rates in MeV^(-1)
table_mu = np.zeros([Enum_init,Enum_fin]) # contains the mu+- -> gamma decay rates dN/dE in MeV^(-1)

table_mu_e = np.zeros([Enum_init,Enum_fin]) # contains the mu+- -> e+- decay rates dN/dE in MeV^(-1)
table_pipm_e = np.zeros([Enum_init,Enum_fin]) # contains the pi+- -> e+- decay rates dN/dE in MeV^(-1)
table_pi0_e = np.zeros([Enum_init,Enum_fin]) # contains the pi0 -> e+- decay rates dN/dE in MeV^(-1); it is actually 0 but needed for formatting reasons


## Decays into photons

for i in range(len(e_init)): # fill in the tables line by line
    table_pi0[i,:] = dnde_pi0(e_fin,e_init[i])
    table_mu[i,:] = dnde_mu(e_fin,e_init[i])
    table_pipm[i,:] = dnde_pi(e_fin,e_init[i])

photon_file = "photon_name.txt" # put in the desired path to your table
print("%15s%15s%15s%15s%15s\n"%("initial_energy","final_energy","pi0>gamma","pi+/->gamma","mu>gamma"),file=open(photon_file,"w"),end = "")
for i in range(Enum_init):
    for j in range(Enum_fin):
        print("%15.5e%15.5e%15.5e%15.5e%15.5e\n"%(e_init[i]/1000.,e_fin[j]/1000.,table_pi0[i,j]*1000.,table_pipm[i,j]*1000.,table_mu[i,j]*1000.),file=open(photon_file,"a"),end = "") # we convert back in GeV and GeV^(-1) for the energies and rates for BlackHawk

## Decays into electrons

for i in range(len(e_init)): # fill in the tables line by line
    table_mu_e[i,:] = dnde_p_mu(e_fin,e_init[i])
    table_pipm_e[i,:] = dnde_p_pi(e_fin,e_init[i])

electron_file = "electron_name.txt" # put in the desired path to your table
print("%15s%15s%15s%15s%15s\n"%("initial_energy","final_energy","pi0>e+/-","pi+/->e+/-","mu+/->e+/-"),file=open(electron_file,"w"),end = "")
for i in range(Enum_init):
    for j in range(Enum_fin):
        print("%15.5e%15.5e%15.5e%15.5e%15.5e\n"%(e_init[i]/1000.,e_fin[j]/1000.,table_pi0_e[i,j]*1000.,table_pipm_e[i,j]*1000.,table_mu_e[i,j]*1000.),file=open(electron_file,"a"),end = "")
