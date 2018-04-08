#include "stdafx.h"
#include "QAS.h"


QAS::QAS()
{

}

QAS::QAS(int n, int m)
{
	N = n;
	M = m;
	spinsInit();
	JInit();

}

QAS::QAS(int n, int m, double t, double gammaa, int nt, int snapnt, double kb )
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

QAS::QAS(int n, int m, double t, double gammaa, int nt, int snapnt, double kb, int gsteps, int tsteps)
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

	Gsteps = gsteps;
	Tsteps = tsteps;
}

void QAS::flipSpin() {
	//choose random spin on NxM grid
	int tmpi, tmpk;
	//srand(time(NULL));
	tmpi = rand() % N;
	tmpk = rand() % M;
	
	//calculates the energy difference for the Metropolis algorithm
	double deltaEnergy = energyDifferenceNN(tmpi, tmpk);
	//double deltaEnergy = energyDifferenceArbitraryJ(tmpi, tmpk);
	
	//Metropolis algorithm
	if (deltaEnergy < 0)
	{
		//cout << "O";
		spins[tmpi][tmpk] *= -1; //* spins[tmpi][tmpk];
		//keeping track of magnetization
		magnetization += 2 * spins[tmpi][tmpk];
	} else 
	{
		double boltzman = exp(- deltaEnergy * beta);
		double uniform = rand() /(double) RAND_MAX;
		if (uniform < boltzman) {
			spins[tmpi][tmpk] *= -1;
			//cout << "_";
			//keeping track of magnetization
			magnetization += 2 * spins[tmpi][tmpk];
		}
	}
}

void QAS::printMagnetization()
{
	cout << "Magnetization = " << magnetization/(M*N) << endl;
}

void QAS::printMagAvgTime()
{
	cout << "MagAvgTime = " << magAvgTime << endl;
}

double QAS::energyDifferenceArbitraryJ(int n, int m)
{
	int tmp = 0;
	for (size_t i = 0; i < N; i++)
	{
		tmp += J[n][i] * spins[i][m];
	}

	if (n != 0 && m != 0)
		return 2.0 *(spins[n][m] * tmp/ ((double)M) + 0.5 / beta * log(1.0 / tanh(beta*gamma / M)) *
			spins[n][m] * (spins[n][m - 1] + spins[n][(m + 1) % M]));
	else if (n == 0 && m != 0)
		return 2.0 * (spins[n][m] * tmp/ ((double)M) + 0.5 / beta * log(1.0 / tanh(beta*gamma / M)) *
			spins[n][m] * (spins[n][m - 1] + spins[n][(m + 1) % M]));
	else if (m == 0 && n != 0)
		return 2.0 * (spins[n][m] * tmp/ ((double)M) + 0.5 / beta * log(1.0 / tanh(beta*gamma / M)) *
			spins[n][m] * (spins[n][M - 1] + spins[n][(m + 1) % M]));
	else if (n == 0 && m == 0)
		return 2.0 * (spins[n][m] * tmp/ ((double)M) + 0.5 / beta * log(1.0 / tanh(beta*gamma / M)) *
			spins[n][m] * (spins[n][M - 1] + spins[n][(m + 1) % M]));
}

double QAS::energyDifferenceNN(int n, int m)
{
	if(n != 0 && m != 0)
		return 2.0 *(spins[n][m] * (J[n][n - 1] * spins[n - 1][m] + J[n][(n + 1) % N] *
			spins[(n + 1) % N][m]) /( (double)M) + 0.5 / beta * log(1.0 / tanh(beta*gamma / M)) *
			spins[n][m] * (spins[n][m - 1] + spins[n][(m + 1) % M]));
	else if(n == 0 && m != 0)
		return 2.0 * (spins[n][m] * (J[n][N - 1] * spins[N - 1][m] + J[n][(n + 1) % N] *
			spins[(n + 1) % N][m]) / ((double)M) + 0.5 / beta * log(1.0 / tanh(beta*gamma / M)) *
			spins[n][m] * (spins[n][m - 1] + spins[n][(m + 1) % M]));
	else if(m == 0 && n != 0)
		return 2.0 * (spins[n][m] * (J[n][n - 1] * spins[n - 1][m] + J[n][(n + 1) % N] *
			spins[(n + 1) % N][m]) / ((double)M) + 0.5 / beta * log(1.0 / tanh(beta*gamma / M)) *
			spins[n][m] * (spins[n][M - 1] + spins[n][(m + 1) % M]));
	else if ( n == 0 && m == 0)
		return 2.0 * (spins[n][m] * (J[n][N - 1] * spins[N - 1][m] + J[n][(n + 1) % N] *
			spins[(n + 1) % N][m]) / ((double)M) + 0.5 / beta * log(1.0 / tanh(beta*gamma / M)) *
			spins[n][m] * (spins[n][M - 1] + spins[n][(m + 1) % M]));
}

void QAS::simulation() {
	spinsInit();
	magHistory.clear();
	for (size_t i = 0; i < NT; i++)
	{
		flipSpin();
		if (i % snapNT == 0) {
			magHistory.push_back(magnetization / (N*M));
		}
	}
	magAvgOverTime();
}

void QAS::simulation(double t) {
	magHistory.clear();
	for (size_t i = 0; i < NT; i++)
	{
		flipSpin();		
		magHistory.push_back(magnetization / (N*M));
	}

}

void QAS::simulationMag(double t) {

	for (size_t i = 0; i < NT; i++)
	{
		flipSpin();
	}

}

double QAS::simulation(double t, double gammaa) {
	T = t;
	gamma = gammaa;
	spinsInit();
	magHistory.clear();
	for (size_t i = 0; i < NT; i++)
	{
		flipSpin();
		if (i % snapNT == 0) {
			magHistory.push_back(magnetization / (N*M));
		}
	}
	magAvgOverTime();
	return magAvgTime;
}

void QAS::spinsInit() {
	srand(time(NULL));
	magnetization = 0;
	this->spins.resize(this->N, vector<int>(this->M, 0));
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < M; j++)
		{
			spins[i][j] = 2 *( rand() % 2 )- 1;
			magnetization += spins[i][j];
		}
	}
}

void QAS::JInit()
{
	this->J.resize(this->N, vector<int>(this->N, 0));
	for (size_t i = 0; i < N - 1; i++)
	{
		J[i][i + 1] = 1;
		J[i + 1][i] = 1;
	}
	J[0][N - 1] = 1;
	J[N - 1][0] = 1;
}

void QAS::printSpins()
{
	cout << endl;
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < M; j++)
		{
			cout << spins[i][j] << "\t";
		}
		cout << endl;
	}
}

void QAS::printJ()
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

QAS::~QAS()
{
}

void QAS::magAvgOverTime()
{
	double integral = 0;
	int lenEnd = magHistory.size();
	int len0 = (int)(magHistory.size()/2);
	for (size_t i = len0; i < lenEnd; i++)
	{
		integral += magHistory[i];
	}
	magAvgTime = integral / (lenEnd - len0);
		
}

void QAS::simulationGammaRange() {
	gammaHistory.clear();
	gammaHistory.resize(2, vector<double>(this->Gsteps, 0));
	
	vector<double> gammaRange;
	gammaRange.resize(Gsteps, 0);
	double tmp = gamma;
	double dgamma = gamma / ((double)Gsteps);
	
	for (size_t i = 0; i < Gsteps; i++)
	{
		gammaRange[i] = tmp;
		tmp -= dgamma;
	
	}
	double magTmp = 0;
	for (size_t i = 0; i < Gsteps; i++)
	{
		magTmp = simulation(T, gammaRange[i]);
		gammaHistory[0][i] = gammaRange[i];
		gammaHistory[1][i] = abs(magTmp);
	}
		
}

void QAS::printGammaHistory()
{
	for (size_t i = 0; i < Gsteps; i++)
	{
		for (size_t j = 0; j < 2; j++)
		{
			cout << gammaHistory[j][i] << "\t";
		}
		cout << endl;
	}
}

void QAS::saveGammaHistory()
{
	ofstream myfile;
	myfile.open("gammaHistory.txt");
	for (size_t i = 0; i < Gsteps; i++)
	{
		myfile << gammaHistory[0][i] << "\t" << gammaHistory[1][i] << endl;
	}
	myfile.close();
}

void QAS::annealing(double tk) {
	vector<double> tRange;
	tRange.resize(Tsteps, 0);
	double tau = pow(tk / T, 1.0 / (double)Tsteps);

	tRange[0] = T;
	for (size_t i = 1; i < Tsteps; i++)
	{
		tRange[i] = tRange[i - 1] * tau;
	}
	
	simulation();
	for (size_t i = 0; i < Tsteps - 1; i++)
	{
		simulation(tRange[i]);
		//cout << tRange[i] << endl;
	}
	simulationMag(tRange[Tsteps - 1]);
	magAvgOverTime();
}

void QAS::annealing(double tk, double gammaa)
{
	gamma = gammaa;
	annealing(tk);
}

void QAS::annealigTempGammaRange(int nt, int ngamma)
{
	annealingMatrix.clear();
	annealingMatrix.resize(nt, vector<double>(ngamma, 0));

	vector<double> gammaRange;
	gammaRange.resize(ngamma, 0);
	double tmp = gamma;
	double dgamma = gamma / ((double)ngamma);

	for (size_t i = 0; i < Gsteps; i++)
	{
		gammaRange[i] = tmp;
		tmp -= dgamma;

	}

	vector<double> tRange;
	tRange.resize(nt, 0);
	tmp = T;
	double dt = T / ((double)nt);

	for (size_t i = 0; i < nt; i++)
	{
		tRange[i] = tmp;
		tmp -= dt;

	}

	for (size_t i = 0; i < ngamma; i++)
	{
		for (size_t j = 0; j < nt; j++)
		{
			annealing(tRange[j], gammaRange[i]);
			annealingMatrix[j][i] = abs(magAvgTime);
		}
		cout << i + 1 << " / " << ngamma << endl;
	}

}

void QAS::saveAnnealingHistory(int nt, int ngamma)
{
	ofstream myfile;
	myfile.open("annealingMatrix.txt");
	for (size_t i = 0; i < ngamma; i++)
	{
		for (size_t j = 0; j < nt; j++)
		{
			myfile << annealingMatrix[j][i] << "\t";
		}
		myfile << endl;
	}
	myfile.close();
}

void QAS::printAnnealingHistory(int nt, int ngamma)
{

	for (size_t i = 0; i < ngamma; i++)
	{
		for (size_t j = 0; j < nt; j++)
		{
			cout << annealingMatrix[j][i] << "\t";
		}
		cout << endl;
	}
}