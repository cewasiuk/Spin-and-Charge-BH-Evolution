// main01.cc is a part of the PYTHIA event generator.
// Copyright (C) 2018 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk.
// It studies the charged multiplicity distribution at the LHC.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

const int nb_fin_part = 6;
const int nb_events = 100000;
const int nb_fin_en = 500;
const int nb_init_en = 250;
const double Emin_init = 5.;
const double Emax_init = 100000.;
const double Emin_fin = 1.e-6;
const double Emax_fin = 100000.;
int WIDTH = 20;
int PRECISION = 10;

int main() {
	char *beam_energy = new char [60];
	Pythia pythia;
	pythia.readString("Beams:idA = 11");
	pythia.readString("Beams:idB = -11");
	pythia.readString("PhaseSpace:pTHatMin = 0.");
	pythia.readString("WeakSingleBoson:ffbar2ffbar(s:gmZ) = on");
	pythia.readString("PartonLevel:ISR = off");
	pythia.readString("23:onMode = off");
	pythia.readString("23:onIfAny = 4"); // decay in c cbar
	pythia.readString("211:mayDecay = on");
	pythia.readString("130:mayDecay = on");
	pythia.readString("321:mayDecay = on");
	pythia.readString("2112:mayDecay = on");
	pythia.readString("13:mayDecay = on");
	int **photon = new int *[nb_init_en];
	int **electron = new int *[nb_init_en];
	int **nu_e = new int *[nb_init_en];
	int **nu_mu = new int *[nb_init_en];
	int **nu_tau = new int *[nb_init_en];
	int **proton = new int *[nb_init_en];
	for(int i = 0;i<nb_init_en;i++){
		photon[i] = new int [nb_fin_en];
		electron[i] = new int [nb_fin_en];
		nu_e[i] = new int [nb_fin_en];
		nu_mu[i] = new int [nb_fin_en];
		nu_tau[i] = new int [nb_fin_en];
		proton[i] = new int [nb_fin_en];
		for(int j = 0;j<nb_fin_en;j++){
			photon[i][j] = 0;
			electron[i][j] = 0;
			nu_e[i][j] = 0;
			nu_mu[i][j] = 0;
			nu_tau[i][j] = 0;
			proton[i][j] = 0;
		}
	}
  	double *final_energies = new double [nb_fin_en];
	for(int i = 0;i<nb_fin_en;i++){
		final_energies[i] = pow(10.,log10(Emin_fin) + (log10(Emax_fin) - log10(Emin_fin))/(nb_fin_en-1)*i);
	}
  	int counter = 0;
	double *beam_energies = new double [nb_init_en];
	for(int i = 0;i<nb_init_en;i++){
		beam_energies[i] = pow(10.,log10(Emin_init) + (log10(Emax_init) - log10(Emin_init))/(nb_init_en-1)*i);
  	}
	for(int k = 0;k<nb_init_en;k++){
		cout << k+1 << "/" << nb_init_en << "\t";
		sprintf(beam_energy,"Beams:eCM = %lf",beam_energies[k]*2);
		pythia.readString(beam_energy);
		pythia.init();
		for (int iEvent = 0; iEvent < nb_events; ++iEvent) {
			if (!pythia.next()){
				continue;
			}
			for(int i = 0;i<pythia.event.size();i++){
				if(pythia.event[i].isFinal() && pythia.event[i].id() == 22){ // photon counter
					if(pythia.event[i].e()>Emin_fin && pythia.event[i].e()<Emax_fin){
						counter = 0;
						while(pythia.event[i].e()>final_energies[counter]){
							counter++;
						}
						photon[k][counter-1]++;
					}
				}
				if(pythia.event[i].isFinal() && (pythia.event[i].id() == 11 || pythia.event[i].id() == -11)){ // electron counter
					if(pythia.event[i].e()>Emin_fin && pythia.event[i].e()<Emax_fin){
						counter = 0;
						while(pythia.event[i].e()>final_energies[counter]){
							counter++;
						}
						electron[k][counter-1]++;
					}
				}
				if(pythia.event[i].isFinal() && (pythia.event[i].id() == 12 || pythia.event[i].id() == -12)){ // electron neutrino counter
					if(pythia.event[i].e()>Emin_fin && pythia.event[i].e()<Emax_fin){
						counter = 0;
						while(pythia.event[i].e()>final_energies[counter]){
							counter++;
						}
						nu_e[k][counter-1]++;
					}
				}
				if(pythia.event[i].isFinal() && (pythia.event[i].id() == 14 || pythia.event[i].id() == -14)){ // muon neutrino counter
					if(pythia.event[i].e()>Emin_fin && pythia.event[i].e()<Emax_fin){
						counter = 0;
						while(pythia.event[i].e()>final_energies[counter]){
							counter++;
						}
						nu_mu[k][counter-1]++;
					}
				}
				if(pythia.event[i].isFinal() && (pythia.event[i].id() == 16 || pythia.event[i].id() == -16)){ // tau neutrino counter
					if(pythia.event[i].e()>Emin_fin && pythia.event[i].e()<Emax_fin){
						counter = 0;
						while(pythia.event[i].e()>final_energies[counter]){
							counter++;
						}
						nu_tau[k][counter-1]++;
					}
				}
				if(pythia.event[i].isFinal() && (pythia.event[i].id() == 2212 || pythia.event[i].id() == -2212)){ // proton counter
					if(pythia.event[i].e()>Emin_fin && pythia.event[i].e()<Emax_fin){
						counter = 0;
						while(pythia.event[i].e()>final_energies[counter]){
							counter++;
						}
						proton[k][counter-1]++;
					}
				}
			}
		}
	}
	ofstream file("table_charm.txt",ios::out);
	file 	<< setw(WIDTH) << "init_energy" << setw(WIDTH) << "fin_energy"
			<< "\t" << setw(WIDTH) << "photon"
			<< "\t" << setw(WIDTH) << "electron"
			<< "\t" << setw(WIDTH) << "nu_e"
			<< "\t" << setw(WIDTH) << "nu_mu"
			<< "\t" << setw(WIDTH) << "nu_tau"
			<< "\t" << setw(WIDTH) << "proton" << endl;
	for(int i = 0;i<nb_init_en;i++){
		for(int j = 0;j<nb_fin_en;j++){
			file << setw(WIDTH) << setprecision(PRECISION) << beam_energies[i] << setw(WIDTH) << setprecision(PRECISION) << final_energies[j];
			file	<< "\t" << setw(WIDTH) << (float)photon[i][j]/nb_events
					<< "\t" << setw(WIDTH) << (float)electron[i][j]/nb_events
					<< "\t" << setw(WIDTH) << (float)nu_e[i][j]/nb_events
					<< "\t" << setw(WIDTH) << (float)nu_mu[i][j]/nb_events
					<< "\t" << setw(WIDTH) << (float)nu_tau[i][j]/nb_events
					<< "\t" << setw(WIDTH) << (float)proton[i][j]/nb_events
					<< endl;
  		}
	}
  	file.close();
	return 0;
}
