#include "stdafx.h"
#include "QTS.h"


QTS::QTS()
{
}

QTS::QTS(int n, int m, double t, double gammaa, int nt, int snapnt, double kb)
{
	N = n;
	M = m;
	T = t;
	gamma = gammaa;
	NT = nt;
	snapNT = snapnt;
	kB = kb;
	beta = 1 / (T * kB);
	spinsInit();
	JInit();
}

QTS::QTS(int n, int m, double t, double gammaa, int nt, int snapnt, int gsteps)
{
	N = n;
	M = m;
	T = t;
	gamma = gammaa;
	Gsteps = gsteps;
	NT = nt;
	snapNT = snapnt;
	beta = 1 / (T * kB);
	spinsInit();
	JInit();
}



QTS::~QTS()
{
}

void QTS::JInit()
{
	srand(time(NULL));
	this->J.resize(this->N, vector<double>(this->N, 0));
	
	locationsX.resize(N, 0);
	locationsY.resize(N, 0);

	for (size_t i = 0; i < N; i++)
	{
		locationsX[i] = rand() / ((double)RAND_MAX);
		locationsY[i] = rand() / ((double)RAND_MAX);
	}
	
	double dij = 0;
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			dij = sqrt(pow((locationsX[j] - locationsX[i]), 2) + pow((locationsY[j] - locationsY[i]), 2));
			J[i][j] = dij;
			J[j][i] = dij;
		}
	}


	
}

void QTS::printJ()
{
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			cout << J[i][j] << "\t";
		}
		cout << endl;
	}
}

void QTS::saveJ()
{
	ofstream myfile;
	myfile.open("J.txt");
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			myfile << J[i][j] << "\t";
		}
		myfile << endl;
	}
	myfile.close();
}

void QTS::spinsInit() {
	srand(time(NULL));
	
	int newsize = N*N;
	this->spins.resize(newsize, vector<int>(this->M, 0));

	vector< vector<int> > spinsHelp;
	spinsHelp.resize(newsize, vector<int>(M, 0));

	int tmp;
	int tmp2;

	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < M; j++)
		{
			do
			{
				tmp2 = rand() % N;
				tmp = N * i + tmp2 ;
				if (spinsHelp[tmp][j] == 0) {
					for (size_t l = 0; l < N; l++)
					{
						spinsHelp[l * N + tmp2][j] = 1;
					}
					spins[tmp][j] = 1;
					break;
				}
			} while (1);
		}
	}
}

void QTS::printSpins()
{
	cout << endl;
	int newsize = N * N;
	for (size_t i = 0; i < M; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			for (size_t k = 0; k < N; k++)
			{
				cout << spins[k*N + j][i] << "\t";

			}
			cout << endl;
		}
		cout <<endl<< "---------------------";
		cout << endl;
	}
}

void QTS::saveSpins()
{
	ofstream myfile;
	myfile.open("pahtsTaken.txt");
	int newsize = N * N;
	for (size_t i = 0; i < M; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			for (size_t k = 0; k < N; k++)
			{
				myfile << spins[k*N + j][i] << "\t";

			}
			myfile << endl;
		}
		//cout << endl;
	}
}

void QTS::simulation() {
	srand(time(NULL));
	for (size_t i = 0; i < NT; i++)
	{
		flipSpin();
	}

}

void QTS::flipSpin() {
	//choose random spin on N*t xM grid
	int tmpi, tmpk;
	
	do
	{
		tmpi = rand() % (N*N);
		tmpk = rand() % M;
	} while (spins[tmpi][tmpk] != 0);
	
	//cout << spins[tmpi][tmpk] << "\t";
	//calculates the energy difference for the Metropolis algorithm
	double deltaEnergy = energyDifference(tmpi, tmpk);
	//cout << deltaEnergy << endl;
	//cout << tmpi << " " << tmpk << endl;

	//Metropolis algorithm
	if (deltaEnergy < 0)
	{
		changeSpins(tmpi, tmpk);
		//changeSpinsBack(tmpi, tmpk);
	}
	else
	{
		double boltzman = exp(-deltaEnergy / T);
		double uniform = rand() / (double)RAND_MAX;
		if (uniform < boltzman) {
			changeSpins(tmpi, tmpk);
		}
	}
}

void QTS::changeSpins(int tmpi, int tmpk)
{
	int t = tmpi % N;
	int n = tmpi / N;
	int where1;
	//cout << "O";
	for (size_t i = 0; i < N; i++)
	{
		if (spins[i * N + t][tmpk] == 1)
		{
			where1 = i;
			spins[i * N + t][tmpk] = 0;
		}
	}

	for (size_t i = 0; i < N; i++)
	{
		if (spins[n * N + i][tmpk] == 1)
		{
			spins[n*N + i][tmpk] = 0;
			spins[N*where1 + i][tmpk] = 1;
			whereIs1 = N * where1 + i;
			break;
		}
	}
	spins[tmpi][tmpk] = 1;

}

void QTS::changeSpinsBack(int tmpi, int tmpk)
{
	spins[whereIs1][tmpk] = 0;
	spins[tmpi][tmpk] = 0;
	spins[(whereIs1 / N) * N + tmpi % N][tmpk] = 1;
	spins[whereIs1 % N + (tmpi / N) * N][tmpk] = 1;
}

double QTS::energyDifference(int tmpi, int tmpk)
{
	double enBefore = spinsEnergy(tmpi, tmpk);
	changeSpins(tmpi, tmpk);
	double enAfter = spinsEnergy(tmpi, tmpk);
	changeSpinsBack(tmpi, tmpk);

	return enAfter - enBefore;
}

double QTS::spinsEnergy(int tmpi, int k)
{
	double energy = 0;
	int iCurr = 0;
	int iNext = 0;
	for (size_t t = 0; t < N; t++)
	{
		for (size_t i = 0; i < N; i++)
		{
			if (spins[N*i + t][k] == 1)
			{
				iCurr = i;
				break;
			}
		}

		for (size_t i = 0; i < N; i++)
		{
			if (spins[N*i +(( t + 1)%N)][k] == 1)
			{
				iNext = i;
				break;
			}
		}
		energy += J[iNext][iCurr];
	}
	if (k > 0)
	{
		energy += 2 * T*log(1 / tanh(gamma / (M*T)))*spins[tmpi][k] * (spins[tmpi][k - 1] + spins[tmpi][(k + 1) % M]);
	}
	else
	{
		energy += 2 * T*log(1 / tanh(gamma / (M*T)))*spins[tmpi][k] * (spins[tmpi][M-1] + spins[tmpi][(k + 1)]);
	}
	//cout << energy;
	return energy;
}

void QTS::printPathLength() {
	double energy = 0;
	int iCurr = 0;
	int iNext = 0;
	double bestPathLen = 10000000;
	for (size_t k = 0; k < M; k++)
	{
		energy = 0;
		iCurr = 0;
		iNext = 0;
		for (size_t t = 0; t < N; t++)
		{
			for (size_t i = 0; i < N; i++)
			{
				if (spins[N*i + t][k] == 1)
				{
					iCurr = i;
					break;
				}
			}

			for (size_t i = 0; i < N; i++)
			{
				if (spins[N*i + ((t + 1) % N)][k] == 1)
				{
					iNext = i;
					break;
				}
			}
			energy += J[iNext][iCurr];
		}
		if (energy < bestPathLen) {
			bestPathLen = energy;
			bestPath = k;
		}
		cout << "Len. of path " << k << ": " << energy << endl;
	}
	cout << endl << "Best path is: " << bestPath << ". Len. of best path: " << bestPathLen << endl;
}

void QTS::simulationGammaRange(double gammak) {
	vector<double> gammaRange;
	gammaRange.resize(Gsteps, 0);
	double tau = pow(gammak / gamma, 1.0 / (double)Gsteps);

	gammaRange[0] = gamma;
	for (size_t i = 1; i < Gsteps; i++)
	{
		gammaRange[i] = gammaRange[i - 1] * tau;
	}
	
	for (size_t i = 0; i < Gsteps; i++)
	{
		simulation();
		cout << i << " / " << Gsteps << endl;
	}
}

void QTS::saveLocations() {
	ofstream myfile;
	myfile.open("points.txt");
	
	for (size_t i = 0; i < locationsX.size(); i++)
	{
		myfile << locationsX[i] << "\t" << locationsY[i] << endl;
	}
	myfile.close();
}

void QTS::saveBestPath() {
	ofstream myfile;
	myfile.open("bestPath.txt");

	for (size_t j = 0; j < N; j++)
	{
		for (size_t k = 0; k < N; k++)
		{
			myfile << spins[k*N + j][bestPath] << "\t";

		}
		myfile << endl;
	}
	myfile.close();
}