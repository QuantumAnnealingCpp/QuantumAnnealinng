#pragma once
#ifndef _QAS_H_ 
#define _QAS_H_
#include<vector>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;


class QAS
{
public:
	vector< vector<int> > spins; // 2D system of spins
	vector< vector<int> > J;
	vector<double> magHistory;
	vector<vector<double>> gammaHistory;
	vector<vector<double>> annealingMatrix;
	double magnetization;
	double T; //thermodynamical temperature
	double gamma; //transverse field strenght factor
	double kB ; //Boltzman constant
	double beta;
	double magAvgTime;
	int N; //N quantum spins
	int M; //additional dimension of spins
	int NT; //liczba kroków czasowych, w których losowany jest jeden spin
	int snapNT; //a value of magnetization is saved every 1000 steps of the simulation
	int Gsteps; //number of evaluation points of gamma for the | <s> | = f(gamma) plot
	int Tsteps; //number of times temperature is decreased during annealing
	void spinsInit(); //initialization of grid of spins
	void JInit(); //initialization of grid of interactions
	void printSpins(); //prints the array of spins
	void printJ(); //prints the array of iteractions
	QAS();
	QAS(int n, int m);
	QAS(int n, int m, double t, double gammaa, int nt, int snapnt, double kb);
	QAS(int n, int m, double t, double gammaa, int nt, int snapnt, double kb, int gsteps, int tsteps);
	void flipSpin(); //decides whether to change the state of a random spin according to Metropolis algorithm
	void printMagnetization(); //prints the magnetization field
	void printMagAvgTime(); //prints the averaged (over time) magnetization
	double energyDifferenceArbitraryJ(int n, int m); //calculates the energy difference for arbitrary interaction matrices J
	void simulation(); //one simulation (all parameters constant)
	void simulation(double t); //different functionality
	void simulationMag(double t); //different functionality
	double simulation(double T, double gamma); //simulation for specified temperature and gamma
	//function uses Monthe Carlo method for searching the minimum of the effective hamiltonian
	double energyDifferenceNN(int n, int m); //calculates the energy difference taking into consideration only the nearest neighbours
	~QAS();
	void magAvgOverTime(); //averages magnetization over time
	void simulationGammaRange(); //simulation with varing gamma
	void printGammaHistory(); //prints the gammaHistory matrix 
	void saveGammaHistory(); //saves into file the gammaHistory matrix
	void annealing(double tk = 0.05); //simulates annealing (varing temperature)
	void annealing(double t, double gammaa);
	void annealigTempGammaRange(int nt, int ngamma);
	void saveAnnealingHistory(int nt, int ngamma);
	void printAnnealingHistory(int nt, int ngamma);
};

#endif 