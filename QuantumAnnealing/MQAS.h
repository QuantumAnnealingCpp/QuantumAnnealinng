#pragma once
#ifndef _MQAS_H_ 
#define _MQAS_H_
#include<vector>
#include<fstream>
#include"QAS.h"
class MQAS
{
public:
	vector<QAS> mqas;
	vector<vector<double>> gammaHistory;
	int numberOfSystems;
	MQAS();
	MQAS(int NumberOfSystems, int M, double T, double gamma, int NT, int snapNT, double kB, int Gsteps, int Tsteps);
	void NRangeGamma();
	~MQAS();
	void printNRangeGamma();
	void saveNRangeGamma();
};

#endif 