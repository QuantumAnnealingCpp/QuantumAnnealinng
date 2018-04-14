
// QuantumAnnealing.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "QAS.h"
#include "MQAS.h"
#include "functions.h"
#include <iostream>
#include "QTS.h"

int main()
{
	double T		= 0.1; //thermodynamical temperature
	double gamma	= 2.0;//transverse field strenght factor
	double kB		= 1.0; //Boltzman constan
	double beta		= 1/kB/T;
	int	N			= 32; //N quantum spins
	int M			= 32; //additional dimension of spins
	int NT			= 2000000; //liczba kroków czasowych, w których losowany jest jeden spin
	int snapNT		= NT/1000; //a value of magnetization is saved every 1000 steps of the simulation
	int Gsteps		= 20; //number of evaluation points of gamma for the | <s> | = f(gamma) plot
	int Tsteps		= 20; //number of times temperature is decreased during annealing
	int numOfSys	= 5; //number of systems of different sizes used to simulate the |magnetization(gamma)| relation
	int nTemperature= 10;
	int nGamma		= 10;

	//QAS qas = QAS(N, M, T, gamma, NT, snapNT, kB,Gsteps, Tsteps); //create instance of Quantum_Annealing System
	
	//Annealing matrix
	/*qas.annealigTempGammaRange(nTemperature, nGamma);
	qas.printAnnealingHistory(nTemperature, nGamma);
	qas.saveAnnealingHistory(nTemperature, nGamma);*/

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


	QTS qts = QTS(N, M, T, gamma, NT, snapNT, Gsteps); //create instance of Quantum_Annealing System
	//qts.printJ();
	//qts.printSpins();
	qts.printPathLength();
	//qts.saveBestPath();
	//qts.simulation();
	qts.simulationGammaRange();
	cout << endl;
	//qts.printSpins();
	qts.printPathLength();
	//qts.spinsEnergy(1, 0);
	qts.saveJ();
	qts.saveSpins();
	qts.saveLocations();
	qts.saveBestPath();

	system("PAUSE");
    return 0;
}

