#pragma once
#ifndef _MQTS_H_ 
#define _MQTS_H_
#include <vector>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include "QTS.h"
using namespace std;

struct params_t
{
	double T;
	double Gamma;
	int M;
	int N;
	int NT;
	int Gsteps;
};

class MQTS
{
public:
	vector<QTS> mqts;
	MQTS();
	MQTS(vector<params_t> &params, int snapnt);

	~MQTS();
	void simulation();
	void simulation(int s);
	void printSpins(int s);
	void saveSpins(int s);
	void saveLocations(int s);
	void saveBestPath(int s);
	void saveJ(int s);
	void MPIsaveJ(int s);
};
#endif