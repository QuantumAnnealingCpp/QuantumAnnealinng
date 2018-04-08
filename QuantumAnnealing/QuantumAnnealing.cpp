
// QuantumAnnealing.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "QAS.h"
#include <iostream>

int main()
{
	double T		= 0.1; //thermodynamical temperature
	double gamma	= 1.5;//transverse field strenght factor
	double kB		= 1.0; //Boltzman constan
	double beta		= 1/kB/T;
	int	N			= 10; //N quantum spins
	int M			= 10; //additional dimension of spins
	int NT			= 100000; //liczba kroków czasowych, w których losowany jest jeden spin
	int snapNT		= NT/1000; //a value of magnetization is saved every 1000 steps of the simulation
	int Gsteps		= 10; //number of evaluation points of gamma for the | <s> | = f(gamma) plot
	int Tsteps		= 30; //number of times temperature is decreased during annealing

	QAS qas = QAS(N, M, T, gamma, NT, snapNT, kB,Gsteps, Tsteps); //create instance of Quantum_Annealing System

	//Check spins and interactions J
	//qas.printJ();
	//qas.printSpins();

	qas.simulation();
	qas.printSpins();
	qas.printMagnetization();
	qas.printMagAvgTime();

	qas.simulationGammaRange();
	qas.printGammaHistory();
	

	system("PAUSE");
    return 0;
}

