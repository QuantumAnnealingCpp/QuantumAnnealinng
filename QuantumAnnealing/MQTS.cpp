#include "stdafx.h"
#include "MQTS.h"
#include <string>
#include <fstream>

MQTS::MQTS()
{
}
MQTS::MQTS(vector<params_t> &params, int snapnt)
{
	mqts.resize(params.size());
	for (int i = 0; i < params.size(); i++)
	{
		mqts[i] = QTS(params[i].N, params[i].M, params[i].T, params[i].Gamma, params[i].NT, snapnt, params[i].Gsteps);
	}
}

void MQTS::simulation()
{
	for (int i = 0; i < mqts.size(); i++)
	{
		mqts[i].printPathLength();
		mqts[i].simulationGammaRange();
		cout << endl;
		saveJ(i);
		saveSpins(i);
		saveLocations(i);
		saveBestPath(i);
	}
}
void MQTS::simulation(int s)
{
	mqts[s].printPathLength();
	mqts[s].simulationGammaRange();
	cout << endl;
	saveJ(s);
	saveSpins(s);
	saveLocations(s);
	saveBestPath(s);
	cout << "done" << endl;
}
MQTS::~MQTS()
{
}
void MQTS::MPIsaveJ(int s)
{
	MPI_File fh;
	string filename = "J" + to_string(s) + ".txt";
	char *cstr = new char[filename.length() + 1];
	strcpy_s(cstr, filename.length() + 1, filename.c_str());
	MPI_File_open(MPI_COMM_SELF, cstr, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
	char bufn[1] = {'\n'};
	for (size_t i = 0; i < mqts[s].N; i++)
	{
		for (size_t j = 0; j < mqts[s].N; j++)
		{
			char buf[42];
			snprintf(buf, 42, "%zd \n", i);
			MPI_File_write(fh, buf, strlen(buf), MPI_CHAR, MPI_STATUS_IGNORE);
		}
		MPI_File_write(fh, bufn, 1, MPI_CHAR, MPI_STATUS_IGNORE);
	}
	MPI_File_close(&fh);
	delete[] cstr;
}
void MQTS::saveJ(int s)
{
	ofstream myfile;
	string filename = "J" + to_string(s) + ".txt";
	cout << filename << endl;
	myfile.open(filename);
	for (size_t i = 0; i < mqts[s].N; i++)
	{
		for (size_t j = 0; j < mqts[s].N; j++)
		{
			myfile << mqts[s].J[i][j] << "\t";
		}
		myfile << endl;
	}
	myfile.close();
}
void MQTS::printSpins(int s)
{
	cout << endl;
	int newsize = mqts[s].N * mqts[s].N;
	for (size_t i = 0; i < mqts[i].M; i++)
	{
		for (size_t j = 0; j < mqts[i].N; j++)
		{
			for (size_t k = 0; k < mqts[i].N; k++)
			{
				cout << mqts[s].spins[k*mqts[i].N + j][i] << "\t";

			}
			cout << endl;
		}
		cout << endl << "---------------------";
		cout << endl; 
	}
}

void MQTS::saveSpins(int s)
{
	ofstream myfile;
	string filename = "pathsTaken" + to_string(s) + ".txt";
	myfile.open(filename);
	cout << filename << endl;

	int newsize = mqts[s].N * mqts[s].N;
	for (size_t i = 0; i < mqts[s].M; i++)
	{
		for (size_t j = 0; j < mqts[s].N; j++)
		{
			for (size_t k = 0; k < mqts[s].N; k++)
			{
				myfile << mqts[s].spins[k*mqts[s].N + j][i] << "\t";

			}
			myfile << endl;
		}
		//cout << endl;
	}
}


void MQTS::saveLocations(int s) {
	ofstream myfile;
	string filename = "points" + to_string(s) + ".txt";
	myfile.open(filename);
	for (size_t i = 0; i < mqts[s].locationsX.size(); i++)
	{
		myfile << mqts[s].locationsX[i] << "\t" << mqts[s].locationsY[i] << endl;
	}
	myfile.close();
}

void MQTS::saveBestPath(int s) {
	ofstream myfile;
	string filename = "bestPath" + to_string(s) + ".txt";
	myfile.open(filename);
	for (size_t j = 0; j < mqts[s].N; j++)
	{
		for (size_t k = 0; k < mqts[s].N; k++)
		{
			myfile << mqts[s].spins[k*mqts[s].N + j][mqts[s].bestPath] << "\t";

		}
		myfile << endl;
	}
	myfile.close();
}