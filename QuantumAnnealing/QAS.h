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
	void printSpins();
	void printJ();
	QAS();
	QAS(int n, int m);
	QAS(int n, int m, double t, double gammaa, int nt, int snapnt, double kb);
	QAS(int n, int m, double t, double gammaa, int nt, int snapnt, double kb, int gsteps, int tsteps);
	void flipSpin(); //decides whether to change the state of a random spin according to Metropolis algorithm
	void printMagnetization();
	void printMagAvgTime();
	double energyDifferenceArbitraryJ(int n, int m);
	void simulation();
	void simulation(double t); //different functionality
	double simulation(double T, double gamma);
	//function uses Monthe Carlo method for searching the minimum of the effective hamiltonian
	double energyDifferenceNN(int n, int m); //calculates the energy difference taking into consideration only the nearest neighbours
	~QAS();
	void magAvgOverTime();
	void simulationGammaRange();
	void printGammaHistory();
	void saveGammaHistory();
	void annealing();
};

#endif 