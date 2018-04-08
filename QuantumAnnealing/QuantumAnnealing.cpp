
// QuantumAnnealing.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "QAS.h"
#include "MQAS.h"
#include "functions.h"
#include <iostream>

int main()
{
	double T		= 1.0; //thermodynamical temperature
	double gamma	= 1.2;//transverse field strenght factor
	double kB		= 1.0; //Boltzman constan
	double beta		= 1/kB/T;
	int	N			= 16; //N quantum spins
	int M			= 16; //additional dimension of spins
	int NT			= 1000000; //liczba kroków czasowych, w których losowany jest jeden spin
	int snapNT		= NT/1000; //a value of magnetization is saved every 1000 steps of the simulation
	int Gsteps		= 10; //number of evaluation points of gamma for the | <s> | = f(gamma) plot
	int Tsteps		= 20; //number of times temperature is decreased during annealing
	int numOfSys	= 5;

	QAS qas = QAS(N, M, T, gamma, NT, snapNT, kB,Gsteps, Tsteps); //create instance of Quantum_Annealing System
	
	//Annealing
	qas.annealing();
	qas.printSpins();
	qas.printMagnetization();

	//Check spins and interactions J
	//qas.printJ();
	//qas.printSpins();

	//Just one simulation
	/*qas.simulation();
	qas.printSpins();
	qas.printMagnetization();
	qas.printMagAvgTime();*/

	//gamma range simulation
	/*qas.simulationGammaRange();
	qas.printGammaHistory();
	qas.saveGammaHistory();*/
	
	cout << endl;

	//gamma range simulation for systems of different size
	/*MQAS mqas = MQAS(numOfSys, M, T, gamma, NT, snapNT, kB, Gsteps, Tsteps);
	mqas.NRangeGamma();
	mqas.printNRangeGamma();
	mqas.saveNRangeGamma();*/

	system("PAUSE");
    return 0;
}

