#pragma once
#ifndef _QTS_H_ 
#define _QTS_H_
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;
class QTS
{
public:
	vector< vector<int> > spins; // 2D system of spins
	vector< vector<double> > J;
	vector<double> locationsX;
	vector<double> locationsY;
	vector<double> magHistory;
	vector<vector<double>> gammaHistory;
	vector<vector<double>> annealingMatrix;
	double magnetization;
	double T; //thermodynamical temperature
	double gamma; //transverse field strenght factor
	double kB; //Boltzman constant
	double beta;
	double magAvgTime;
	int N; //N quantum "cities"
	int M; //additional Suzuki Trotter dimension
	int NT; //liczba kroków czasowych, w których losowany jest "miasto"
	int snapNT; //a value of magnetization is saved every 1000 steps of the simulation
	int Gsteps; //number of evaluation points of gamma for the | <s> | = f(gamma) plot
	int Tsteps; //number of times temperature is decreased during annealing
	int whereIs1;

	QTS();
	QTS(int n, int m, double t, double gammaa, int nt, int snapnt, double kb);
	QTS(int n, int m, double t, double gammaa, int nt, int snapnt, int gsteps);
	void spinsInit();
	void printSpins();
	void simulation();
	void flipSpin();
	void changeSpins(int tmpi, int tmpk);
	void changeSpinsBack(int tmpi, int tmpk);
	double energyDifference(int tmpi, int tmpk);
	double spinsEnergy(int tmpi, int k);
	void printPathLength();
	void simulationGammaRange(double gammak = 0.05);
	~QTS();
	void printJ();
	void JInit();
};

#endif
