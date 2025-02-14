import numpy as np
import pylab

Isatis_path = "your path" # here put your path to the Isatis folder
results_name = "your name" # here put the name of your Isatis results

## Defining the plot format

fig_width_pt = 85*2.83465  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height_pt = 80*2.83465
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


## Plotting all constraints
constraints_file = np.genfromtxt("%s/results_photons_%s.txt"%(Isatis_path,results_name),dtype = "str")
constraints_names_bis = constraints_file[0,1:]
constraints = np.zeros([len(constraints_file)-1,len(constraints_file[0])-1])
for i in range(len(constraints)):
    for j in range(len(constraints[0])):
        constraints[i,j] = float(constraints_file[i+1,j+1])

Mmin = # minimum PBH mass in your runs
Mmax = # maximum PBH mass in your runs
masses = np.logspace(np.log10(Mmin),np.log10(Mmax),len(constraints))

# creating labels
constraints_names = []
for i in range(len(constraints_names_bis)):
    temp = constraints_names_bis[i].split("_")
    temp2 = ""
    for i in range(len(temp)-1):
        temp2 = "".join([temp2,temp[i],'\,\,'])
    temp2 = "".join([temp2,'\,\,[arXiv:',temp[-1],']'])
    constraints_names.append(temp2)

f = pylab.figure(1)
f.clf()
ax = f.add_subplot(111)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("$M_{\\rm PBH}\,\,{\\rm (g)}$")
ax.set_ylabel("$f_{\\rm PBH}$")
ax.set_xlim(Mmin,Mmax)
ax.set_ylim(np.min(constraints),2.)
for i in range(len(constraints_names)):
    ax.plot(masses,constraints[:,i],label = "${\\rm %s}$"%(constraints_names[i]))
ax.legend(loc = "best")
ax.grid(True)
f.tight_layout(rect = [0,0,1,1])
f.show()
f.savefig("%s/results.pdf"%(Isatis_path))

## IGRB (all)
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

HEAO1 = np.genfromtxt("%s/constraints/photons/flux_HEAO_balloon_9903492.txt"%(Isatis_path),skip_header = 6)
COMPTEL = np.genfromtxt("%s/constraints/photons/flux_COMPTEL_1502.06116.txt"%(Isatis_path),skip_header = 6)
EGRET = np.genfromtxt("%s/constraints/photons/flux_EGRET_0405441.txt"%(Isatis_path),skip_header = 6)
Fermilat = np.genfromtxt("%s/constraints/photons/flux_Fermi-LAT_1410.3696.txt"%(Isatis_path),skip_header = 6)

energies1 = np.logspace(np.log10(6e-5),2,100)
energies1_bis = np.logspace(np.log10(3e-6),np.log10(6e-5),100)
background1 = (0.0259*(energies1*1e+6/60)**(-5.5) + 0.504*(energies1*1e+6/60)**(-1.58) + 0.0288*(energies1*1e+6/60)**(-1.05))/energies1
background1_bis = 7.877*(energies1_bis*1e+6)**(-0.29)*np.exp(-(energies1_bis*1e+6/41.13))/energies1_bis

energies2 = np.logspace(np.log10(8e-4),1,100)
background2 = 2.74e-6*energies2**(-2)

energies3 = np.logspace(np.log10(2e-5),np.log10(5e-3),100)
background3 = 64.2e+3/((energies3/35.6966e-6)**(1.4199) + (energies3/35.6966e-6)**(2.8956))

energies4 = np.logspace(np.log10(150e-6),np.log10(5e-3),100)
background4 = 0.0044135e+3*(energies4/1e-3)**(-2.8956)

energies5 = np.logspace(np.log10(1.5e-6),np.log10(200e-6),100)
background5 = 0.109e+6/((energies5/(29.0e-6))**(1.40) + (energies5/(29.0e-6))**(2.88))

f = pylab.figure(1)
f.clf()
ax = f.add_subplot(111)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("$E\,\,{\\rm (GeV)}$")
ax.set_ylabel("${\\rm d}\Phi_{\\rm egal}/{\\rm d}E\,\,{\\rm (GeV}^{-1}\cdot{\\rm cm}^{-2}\cdot{\\rm s}^{-1}\cdot{\\rm sr}^{-1})$")
ax.set_ylim(1e-17,1e+7)
ax.set_xlim(1e-6,1e+4)
ax.scatter(HEAO1[:,0],HEAO1[:,3],marker = "+",s = 30,label = "${\\rm HEAO1}$")
ax.scatter(COMPTEL[:,0],COMPTEL[:,3],marker = "+",s = 40,label = "${\\rm COMPTEL}$")
ax.scatter(EGRET[:,0],EGRET[:,3],marker = "+",s = 30,label = "${\\rm EGRET}$")
ax.scatter(Fermilat[:,0],Fermilat[:,3],marker = "+",s = 40,label = "${\\rm Fermi-LAT}$")
ax.plot(energies1_bis,background1_bis,color = "black",linewidth = 1,linestyle = "--")
ax.plot(energies1,background1,color = "black",linewidth = 1,linestyle = "--")
ax.plot(energies2,background2,color = "black",linewidth = 1,linestyle = "-.")
ax.plot(energies3,background3,color = "black",linewidth = 1,linestyle = ":")
ax.plot(energies4,background4,color = "black",linewidth = 1,linestyle = "-")
ax.plot(energies5,background5,color = "black",linewidth = 1,linestyle = "-.")
ax.grid(True)
ax.legend(loc = "best")
f.tight_layout(rect = [0,0,1,1])
f.show()
f.savefig("%s/IGRB.pdf"%(Isatis_path))

f = pylab.figure(2)
f.clf()
ax = f.add_subplot(111)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("$E\,\,{\\rm (GeV)}$")
ax.set_ylabel("${\\rm d}\Phi_{\\rm egal}/{\\rm d}E\,\,{\\rm (GeV}^{-1}\cdot{\\rm cm}^{-2}\cdot{\\rm s}^{-1}\cdot{\\rm sr}^{-1})$")
ax.set_ylim(1e-17,1e+7)
ax.set_xlim(1e-6,1e+4)
ax.errorbar(HEAO1[:,0],HEAO1[:,3],xerr = [HEAO1[:,1],HEAO1[:,2]],yerr = [HEAO1[:,4],HEAO1[:,5]],fmt = "o",markersize = 2,fillstyle = None,label = "${\\rm HEAO1}$")
ax.errorbar(COMPTEL[:,0],COMPTEL[:,3],xerr = [COMPTEL[:,1],COMPTEL[:,2]],yerr = [COMPTEL[:,4],COMPTEL[:,5]],fmt = "o",markersize = 2,fillstyle = None,label = "${\\rm COMPTEL}$")
ax.errorbar(EGRET[:,0],EGRET[:,3],xerr = [EGRET[:,1],EGRET[:,2]],yerr = [EGRET[:,4],EGRET[:,5]],fmt = "o",markersize = 2,fillstyle = None,label = "${\\rm EGRET}$")
ax.errorbar(Fermilat[:,0],Fermilat[:,3],xerr = [Fermilat[:,1],Fermilat[:,2]],yerr = [Fermilat[:,4],Fermilat[:,5]],fmt = "o",markersize = 2,fillstyle = None,label = "${\\rm Fermi-LAT}$")
ax.plot(energies1_bis,background1_bis,color = "black",linewidth = 1,linestyle = "--")
ax.plot(energies1,background1,color = "black",linewidth = 1,linestyle = "--")
ax.plot(energies2,background2,color = "black",linewidth = 1,linestyle = "-.")
ax.plot(energies3,background3,color = "black",linewidth = 1,linestyle = ":")
ax.plot(energies4,background4,color = "black",linewidth = 1,linestyle = "-")
ax.plot(energies5,background5,color = "black",linewidth = 1,linestyle = "-.")
ax.grid(True)
ax.legend(loc = "best")
f.tight_layout(rect = [0,0,1,1])
f.show()
f.savefig("%s/IGRB_errorbars.pdf"%(Isatis_path))

f = pylab.figure(3)
f.clf()
ax = f.add_subplot(111)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("$E\,\,{\\rm (GeV)}$")
ax.set_ylabel("${\\rm d}\Phi_{\\rm egal}/{\\rm d}E\,\,{\\rm (GeV}^{-1}\cdot{\\rm cm}^{-2}\cdot{\\rm s}^{-1}\cdot{\\rm sr}^{-1})$")
#ax.set_ylim(1e-17,1e+7)
ax.set_xlim(1e-6,1e+4)
ax.errorbar(HEAO1[:,0],HEAO1[:,3]*HEAO1[:,0]**2,xerr = [HEAO1[:,1],HEAO1[:,2]],yerr = [HEAO1[:,4]*HEAO1[:,0]**2,HEAO1[:,5]*HEAO1[:,0]**2],fmt = "o",markersize = 2,fillstyle = None,label = "${\\rm HEAO1}$")
ax.errorbar(COMPTEL[:,0],COMPTEL[:,3]*COMPTEL[:,0]**2,xerr = [COMPTEL[:,1],COMPTEL[:,2]],yerr = [COMPTEL[:,4]*COMPTEL[:,0]**2,COMPTEL[:,5]*COMPTEL[:,0]**2],fmt = "o",markersize = 2,fillstyle = None,label = "${\\rm COMPTEL}$")
ax.errorbar(EGRET[:,0],EGRET[:,3]*EGRET[:,0]**2,xerr = [EGRET[:,1],EGRET[:,2]],yerr = [EGRET[:,4]*EGRET[:,0]**2,EGRET[:,5]*EGRET[:,0]**2],fmt = "o",markersize = 2,fillstyle = None,label = "${\\rm EGRET}$")
ax.errorbar(Fermilat[:,0],Fermilat[:,3]*Fermilat[:,0]**2,xerr = [Fermilat[:,1],Fermilat[:,2]],yerr = [Fermilat[:,4]*Fermilat[:,0]**2,Fermilat[:,5]*Fermilat[:,0]**2],fmt = "o",markersize = 2,fillstyle = None,label = "${\\rm Fermi-LAT}$")
ax.plot(energies1_bis,background1_bis*energies1_bis**2,color = "black",linewidth = 1,linestyle = "--")
ax.plot(energies1,background1*energies1**2,color = "black",linewidth = 1,linestyle = "--")
ax.plot(energies2,background2*energies2**2,color = "black",linewidth = 1,linestyle = "-.")
ax.plot(energies3,background3*energies3**2,color = "black",linewidth = 1,linestyle = ":")
ax.plot(energies4,background4*energies4**2,color = "black",linewidth = 1,linestyle = "-")
ax.plot(energies5,background5*energies5**2,color = "black",linewidth = 1,linestyle = "-.")
ax.grid(True)
ax.legend(loc = "best")
f.tight_layout(rect = [0,0,1,1])
f.show()
f.savefig("%s/IGRB_rescaled.pdf"%(Isatis_path))

## Galactic background (all)
INTEGRAL = np.genfromtxt("%s/constraints/photons/flux_INTEGRAL_1107.0200.txt"%(Isatis_path),skip_header = 6)
COMPTEL = np.genfromtxt("%s/constraints/photons/flux_COMPTEL_1107.0200.txt"%(Isatis_path),skip_header = 6)
EGRET = np.genfromtxt("%s/constraints/photons/flux_EGRET_9811211.txt"%(Isatis_path),skip_header = 6)
FermiLAT = np.genfromtxt("%s/constraints/photons/flux_Fermi-LAT_1101.1381.txt"%(Isatis_path),skip_header = 6)

energies = np.logspace(-4,np.log10(20e-3),100)
BARTELS = np.zeros(len(energies)) # see 1703.02546
for i in range(len(energies)):
    BARTELS[i] = 13.*(energies[i]/1e-3)**(-1.8)*np.exp(-energies[i]/(20e-3))

f = pylab.figure(1)
f.clf()
ax = f.add_subplot(111)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(1e-5,1e+2)
ax.set_ylim(1e-10,1e+5)
ax.scatter(INTEGRAL[:,0],INTEGRAL[:,3],marker = "+",s = 30,label = "${\\rm INTEGRAL}$")
ax.scatter(COMPTEL[:,0],COMPTEL[:,3],marker = "+",s = 40,label = "${\\rm COMPTEL}$")
ax.scatter(EGRET[:,0],EGRET[:,3],marker = "+",s = 40,label = "${\\rm EGRET}$")
ax.scatter(FermiLAT[:,0],FermiLAT[:,3],marker = "+",s = 30,label = "${\\rm Fermi-LAT}$")
ax.plot(energies,BARTELS,color = "black",linewidth = 1,linestyle = "--")
ax.grid(True)
ax.legend(loc = "best")
f.show()
f.tight_layout(rect = [0,0,1,1])
f.savefig("%s/galactic.pdf"%(Isatis_path))

f = pylab.figure(2)
f.clf()
ax = f.add_subplot(111)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(1e-5,1e+2)
ax.set_ylim(1e-10,1e+5)
ax.errorbar(INTEGRAL[:,0],INTEGRAL[:,3],xerr = [INTEGRAL[:,1],INTEGRAL[:,2]],yerr = [INTEGRAL[:,4],INTEGRAL[:,5]],fmt = "o",markersize = 2,fillstyle = None,label = "${\\rm INTEGRAL}$")
ax.errorbar(COMPTEL[:,0],COMPTEL[:,3],xerr = [COMPTEL[:,1],COMPTEL[:,2]],yerr = [COMPTEL[:,4],COMPTEL[:,5]],fmt = "o",markersize = 2,fillstyle = None,label = "${\\rm COMPTEL}$")
ax.errorbar(EGRET[:,0],EGRET[:,3],xerr = [EGRET[:,1],EGRET[:,2]],yerr = [EGRET[:,4],EGRET[:,5]],fmt = "o",markersize = 2,fillstyle = None,label = "${\\rm EGRET}$")
ax.errorbar(FermiLAT[:,0],FermiLAT[:,3],xerr = [FermiLAT[:,1],FermiLAT[:,2]],yerr = [FermiLAT[:,4],FermiLAT[:,5]],fmt = "o",markersize = 2,fillstyle = None,label = "${\\rm Fermi-LAT}$")
ax.plot(energies,BARTELS,color = "black",linewidth = 1,linestyle = "--")
ax.grid(True)
ax.legend(loc = "best")
f.show()
f.tight_layout(rect = [0,0,1,1])
f.savefig("%s/galactic_errorbar.pdf"%(Isatis_path))

f = pylab.figure(3)
f.clf()
ax = f.add_subplot(111)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(1e-5,1e+2)
ax.set_ylim(1e-6,1e-4)
ax.errorbar(INTEGRAL[:,0],INTEGRAL[:,3]*INTEGRAL[:,0]**2,xerr = [INTEGRAL[:,1],INTEGRAL[:,2]],yerr = [INTEGRAL[:,4]*INTEGRAL[:,0]**2,INTEGRAL[:,5]*INTEGRAL[:,0]**2],fmt = "o",markersize = 2,fillstyle = None,label = "${\\rm INTEGRAL}$")
ax.errorbar(COMPTEL[:,0],COMPTEL[:,3]*COMPTEL[:,0]**2,xerr = [COMPTEL[:,1],COMPTEL[:,2]],yerr = [COMPTEL[:,4]*COMPTEL[:,0]**2,COMPTEL[:,5]*COMPTEL[:,0]**2],fmt = "o",markersize = 2,fillstyle = None,label = "${\\rm COMPTEL}$")
ax.errorbar(EGRET[:,0],EGRET[:,3]*EGRET[:,0]**2,xerr = [EGRET[:,1],EGRET[:,2]],yerr = [EGRET[:,4]*EGRET[:,0]**2,EGRET[:,5]*EGRET[:,0]**2],fmt = "o",markersize = 2,fillstyle = None,label = "${\\rm EGRET}$")
ax.errorbar(FermiLAT[:,0],FermiLAT[:,3]*FermiLAT[:,0]**2,xerr = [FermiLAT[:,1],FermiLAT[:,2]],yerr = [FermiLAT[:,4]*FermiLAT[:,0]**2,FermiLAT[:,5]*FermiLAT[:,0]**2],fmt = "o",markersize = 2,fillstyle = None,label = "${\\rm Fermi-LAT}$")
ax.plot(energies,BARTELS*energies**2,color = "black",linewidth = 1,linestyle = "--")
ax.grid(True)
ax.legend(loc = "best")
f.show()
f.tight_layout(rect = [0,0,1,1])
f.savefig("%s/galactic_rescaled.pdf"%(Isatis_path))




