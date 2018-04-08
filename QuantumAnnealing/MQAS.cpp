#include "stdafx.h"
#include "MQAS.h"


MQAS::MQAS()
{
}

MQAS::MQAS(int numberofsystems, int M, double T, double gamma, int NT, int snapNT, double kB, int Gsteps, int Tsteps)
{
	numberOfSystems = numberofsystems;
	int n0 = 8;
	mqas.clear();
	mqas.resize(numberOfSystems);
	for (size_t i = 0; i < numberOfSystems; i++)
	{
		mqas[i] = QAS(n0,  M, T, gamma, NT, snapNT, kB, Gsteps, Tsteps);
		n0 *= 2;
	}

}
void MQAS::NRangeGamma()
{
	gammaHistory.clear();
	gammaHistory.resize(numberOfSystems + 1, vector<double>(this->mqas[0].Gsteps, 0));
	for (size_t i = 0; i < numberOfSystems; i++)
	{
		mqas[i].simulationGammaRange();
		for (size_t j = 0; j < mqas[i].Gsteps; j++)
		{
			gammaHistory[0][j] = mqas[i].gammaHistory[0][j];
			gammaHistory[i][j] = mqas[i].gammaHistory[1][j];
		}
	}
}

MQAS::~MQAS()
{
}

void MQAS::printNRangeGamma()
{
	for (size_t i = 0; i < mqas[0].Gsteps; i++)
	{
		for (size_t j = 0; j < numberOfSystems; j++)
		{
			cout << gammaHistory[j][i] << "\t";
		}
		cout << endl;
	}
}

void MQAS::saveNRangeGamma()
{
	ofstream myfile;
	myfile.open("NgammaHistory.txt");
	for (size_t i = 0; i < mqas[0].Gsteps; i++)
	{
		for (size_t j = 0; j < numberOfSystems; j++)
		{
			myfile << gammaHistory[j][i] << "\t";
		}
		myfile << endl;
	}
	myfile.close();
}