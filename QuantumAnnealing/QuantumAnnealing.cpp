
// QuantumAnnealing.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "QAS.h"
#include "MQAS.h"
#include "functions.h"
#include <iostream>

int main()
{
	double T		= 0.5; //thermodynamical temperature
	double gamma	= 1.3;//transverse field strenght factor
	double kB		= 1.0; //Boltzman constan
	double beta		= 1/kB/T;
	int	N			= 10; //N quantum spins
	int M			= 10; //additional dimension of spins
	int NT			= 100000; //liczba kroków czasowych, w których losowany jest jeden spin
	int snapNT		= NT/1000; //a value of magnetization is saved every 1000 steps of the simulation
	int Gsteps		= 10; //number of evaluation points of gamma for the | <s> | = f(gamma) plot
	int Tsteps		= 20; //number of times temperature is decreased during annealing
	int numOfSys	= 5; //number of systems of different sizes used to simulate the |magnetization(gamma)| relation
	int nTemperature= 10;
	int nGamma		= 10;

	QAS qas = QAS(N, M, T, gamma, NT, snapNT, kB,Gsteps, Tsteps); //create instance of Quantum_Annealing System
	
	//Annealing matrix
	qas.annealigTempGammaRange(nTemperature, nGamma);
	qas.printAnnealingHistory(nTemperature, nGamma);
	qas.saveAnnealingHistory(nTemperature, nGamma);

	//Annealing
	/*qas.annealing();
	qas.printSpins();
	qas.printMagnetization();*/

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

