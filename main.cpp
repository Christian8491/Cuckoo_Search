/* This code implements "Cuckoo Search via Lévy Flights" in C++11
*
* Params:
* n: number of nests
* Pa: probability that the guest will discover the foreign egg

* numberDomain: dimension of the vector (for this purpose {2, 8})
**/

#include <fstream>
#include <iostream>
#include <vector>
#include <stdlib.h>								// rand()
#include <time.h>								// time()
#include <math.h>								// tgamma()
#include <chrono>
#include <random>								// normal_distribution()
#include <algorithm>							// random_shuffle()

using namespace std;

# define M_PI				3.14159265358979323846
# define N_ITERACIONES				100

/* ============================================ Objective Functions ================================================ */

/* A d-dimensional Test function 
 * @param x: a d-dimensional vector 
 * @return the total sum of each (x_i - 1)^2 */
double testFuncion(vector<double>& x)
{
	double totalSum = 0.0;
	for (size_t i = 0; i < x.size(); ++i) {
		totalSum += pow(x[i] - 1.0, 2);
	}

	return totalSum;
}

/* This function is characterized by forming waves of local maximums and minimums 
 * @param x: a d-dimensional vector 
 * @return the function3 evaluated in this x vector */
double function3(vector<double>& x){
	double summation = 0.0;

	for (size_t i = 0; i < x.size(); ++i) {
		summation += pow(x[i], 2.0);
	}

	double numerator = pow(sin(sqrt(summation)), 2.0) - 0.5;
	double denominator = pow(1.0 + 0.001 * summation, 2.0);

	return 0.5 - numerator / denominator;
}

/* This function is very complex with many local minimums
* @param x: a d-dimensional vector
* @return the Schwefel evaluated in this x vector */
double Schwefel(vector<double>& x){
	int d = x.size();

	double sumatoria1 = 0.0;
	for (size_t i = 0; i < x.size(); ++i) {
		sumatoria1 += x[i] * sin(sqrt(abs(x[i])));
	}

	return 418.9829 * d - sumatoria1;
}

/* This function is characterized by an almost flat region and a large hole in the center.
* @param x: a d-dimensional vector
* @return the Ackley evaluated in this x vector */
double Ackley(vector<double>& x, double a = 20.0, double b = 0.2, double c = 2 * M_PI){
	int d = x.size();

	// First Summation
	double firstSummation = 0.0;
	for (size_t i = 0; i < x.size(); ++i) {
		firstSummation += pow(x[i], 2);
	}
	firstSummation = sqrt(firstSummation / d);

	// Second Summation
	double secondSummation = 0.0;
	for (size_t i = 0; i < x.size(); ++i) {
		secondSummation += cos(c * x[i]);
	}
	secondSummation /= d;

	double totalValue = (-1 * a * exp(-1 * b * firstSummation)) - exp(secondSummation) + a + exp(1);
	return totalValue;
}


/*================================================ Util Functions ==================================================== */

/* Save in a *.txt format the best values */
void saveData(vector<double>& bestValues) 
{
	ofstream myFile;
	myFile.open("data.txt");

	if (myFile.fail()) {
		cerr << "The file was not found" << endl;
		return;
	}

	for (size_t i = 0; i < bestValues.size(); i++) {
		myFile << bestValues[i];
		myFile << "\n";
	}
}

/* Params are initialized as indicated in the paper */
void initializeParams(int& n, double& Pa)
{
	n = 25;
	Pa = 0.25;
}

/* Only generate random values between two numbers */
double fRand(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

/* Function to calculate normal distribution */
void normalDistribution(vector<double>& normalDistr, double& sigma, bool multiplicate)
{
	// construct a trivial random generator engine from a time-based seed:
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
	normal_distribution<double> distribution(0.0, 1.0);

	if (multiplicate == true) {
		for (size_t i = 0; i < normalDistr.size(); ++i) {
			normalDistr[i] = distribution(generator) * sigma;
		}
	}
	else {
		for (size_t i = 0; i < normalDistr.size(); ++i) {
			normalDistr[i] = distribution(generator);
		}
	}

}

/* Auxiliar algorithm for Mantegna */
void mantegna(vector<double>& step, vector<double>& u, vector<double>& v, double& beta)
{
	for (size_t i = 0; i < step.size(); ++i) {
		step[i] = u[i] / pow(abs(v[i]), 1 / beta);
	}
}


/* Application of simple constraints */
vector<double> simpleBounds(vector<double>& s, vector<double>& lowerBound, vector<double>& upperBound)
{
	// Apply the lower bound
	vector<double> nsTemp = s;
	for (size_t i = 0; i < s.size(); ++i) {
		if (nsTemp[i] < lowerBound[i]) nsTemp[i] = lowerBound[i];
	}

	// Apply the upper bound
	for (size_t i = 0; i < s.size(); ++i) {
		if (nsTemp[i] > upperBound[i]) nsTemp[i] = upperBound[i];
	}

	// Update this new move 
	return nsTemp;
}

/* Get cuckoos by ramdom walk - Levy implementation */
vector<vector<double>> getCuckoos(vector<vector<double>>& nest, vector<double>& bestNest, vector<double>& lowerBound, vector<double>& upperBound)
{
	int n = nest.size();
	vector<vector<double>> newNest(n, vector<double>(nest[0].size()));

	// Levy exponent and coefficient
	double beta = 1.5;
	double sigma = (tgamma(1 + beta)* sin(M_PI * beta / 2) / (tgamma((1 + beta) / 2) * beta * pow(2.0, (beta - 1) / 2)));
	double sigmaPow = pow(sigma, 1 / beta);

	for (int i = 0; i < n; ++i) {
		vector<double> s = nest[i];

		// Levy flights by Mantegna's algorithm
		vector<double> u(s.size()), v(s.size()), step(s.size());
		normalDistribution(u, sigmaPow, true);
		normalDistribution(v, sigmaPow, false);
		mantegna(step, u, v, beta);													// fill step variable

		vector<double> stepsize(s.size());
		for (size_t j = 0; j < stepsize.size(); ++j) {
			stepsize[j] = 0.01 * step[j] * (s[j] - bestNest[j]);
		}

		// Actual random walks or flights
		vector<double> auxNormalD(s.size());
		normalDistribution(auxNormalD, sigmaPow, false);

		for (size_t j = 0; j < s.size(); ++j) {
			s[j] = s[j] + stepsize[j] * auxNormalD[j];
		}

		// Apply simple bounds/limits
		nest[i] = simpleBounds(s, lowerBound, upperBound);			// newNest[i] = simpleBounds(s, lowerBound, upperBound);	
	}

	newNest = nest;

	return newNest;
}

/* Find the current best nest and fill bestNest param */
double getBestNest(vector<double>& bestNest, vector<vector<double>>& nest, vector<vector<double>>& newNest, vector<vector<double>>& fitness)
{
	// Evaluating all new solutions
	for (size_t i = 0; i < nest.size(); ++i) {
		double fNew = Ackley(newNest[i]);							// Here goes the objective function
		if (fNew <= fitness[i][0]) {
			fitness[i][0] = fNew;
			nest[i] = newNest[i];
		}
	}

	// Find the current best
	double fitnessMin = fitness[0][0];
	int k = 0;
	for (size_t i = 1; i < fitness.size(); ++i) {
		if (fitness[i][0] < fitnessMin) {
			fitnessMin = fitness[i][0];
			k = i;
		}
	}

	bestNest = nest[k];

	return fitnessMin;
}

/* rearrange the matrix from rows */
vector<vector<double>> randPerm(vector<vector<double>> nest)
{
	int n = nest.size();
	vector<vector<double>> nestCopy = nest;
	vector<int> index(n);
	for (int i = 0; i < n; ++i) {
		index[i] = i;
	}

	// Random permutation the order
	for (int i = 0; i < n; ++i) {
		int j = rand() % (n - i) + i;
		// Swap i and j
		int t = index[j];
		index[j] = index[i];
		index[i] = t;
	}

	for (int i = 0; i < n; i++) {
		nestCopy[i] = nest[index[i]];
	}

	return nestCopy;
}

/* Replace some nests by constructing new solutions/nests */
vector<vector<double>> emptyNests(vector<vector<double>>& nest, vector<double>& lowerBound, vector<double>& upperBound, double& Pa)
{
	int n = nest.size();

	// Discovered or not -- a status vector (K contains 0 or 1)
	unsigned seed2 = std::chrono::system_clock::now().time_since_epoch().count();
	srand(seed2);
	vector<vector<int>> K(n, vector<int>(nest[0].size()));
	for (int i = 0; i < n; ++i) {
		for (size_t j = 0; j < nest[0].size(); ++j) {
			if (fRand(-1.0, 1.0) > Pa) K[i][j] = 1.0;
			else  K[i][j] = 0.0;
		}
	}

	// New solution by biased/selective random walks
	unsigned seed3 = std::chrono::system_clock::now().time_since_epoch().count();
	srand(seed3);
	vector<vector<double>> nestRand1 = randPerm(nest);
	vector<vector<double>> nestRand2 = randPerm(nest);
	vector<vector<double>> stepsize(n, vector<double>(nest[0].size()));

	for (int i = 0; i < n; ++i) {
		for (size_t j = 0; j < nest[0].size(); ++j) {
			stepsize[i][j] = fRand(-1.0, 1.0) * (nestRand1[i][j] - nestRand2[i][j]);
		}
	}

	vector<vector<double>> newNest(n, vector<double>(nest[0].size()));
	// Fill newNest
	for (int i = 0; i < n; ++i) {
		for (size_t j = 0; j < nest[0].size(); ++j) {
			newNest[i][j] = nest[i][j] + stepsize[i][j] * K[i][j];
		}

	}

	for (int i = 0; i < n; ++i) {
		vector<double> s = newNest[i];
		newNest[i] = simpleBounds(s, lowerBound, upperBound);
	}

	return newNest;
}

/* Principal function to Cuckoo Search */
void cuckooSearch(int& n, double& Pa)
{
	double tolerance = 1.0e-5;
	int numberDomain = 2;											// x = (x_1, x_2, ... , x_n)
	int nIterations = 0;											// number of iterations

	// Bounds
	vector<double> lowerBound(numberDomain, -5.0);					// limits of domain
	vector<double> upperBound(numberDomain, 5.0);					// limits of domain

	// Random initial solutions
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	srand(seed);
	vector<vector<double>> nest(n, vector<double>(numberDomain));
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < numberDomain; ++j) {
			nest[i][j] = lowerBound[j] + (upperBound[j] - lowerBound[j]) * fRand(-1.0, 1.0);
		}
	}

	// Get the current best
	vector<vector<double>> fitness(n, vector<double>(1));
	for (size_t i = 0; i < fitness.size(); ++i) {
		fitness[i][0] = 1.0e11;
	}

	vector<double> bestNest;
	double fMin = getBestNest(bestNest, nest, nest, fitness);

	// To save data
	vector<double> bestValues(N_ITERACIONES);

	while (fMin > tolerance && nIterations < N_ITERACIONES) {

		// Generate new solutions (but keep the current best)
		vector<vector<double>> newNest = getCuckoos(nest, bestNest, lowerBound, upperBound);
		vector<double> best;
		double fNew = getBestNest(best, nest, newNest, fitness);

		// Discovery and randomization
		newNest = emptyNests(nest, lowerBound, upperBound, Pa);

		// Evaluate this set of solutions
		fNew = getBestNest(best, nest, newNest, fitness);
		cout << "fmin N " << nIterations << ": " << fNew << endl;

		// Find the best objective so far
		if (fNew < fMin) {
			fMin = fNew;
			bestNest = best;
		}
		
		// Update the counter
		nIterations++;
	}

	cout << "\nTotal iterations: " << nIterations << endl;

	cout << "\nminimum value in function: " << fMin << endl;

	cout << "\nBest Nest vector: " << endl;
	for (size_t i = 0; i < bestNest.size(); ++i) {
		cout << bestNest[i] << " ";
	}
	cout << endl;
}

int main()
{
	// Initial parameters
	int n;							// number of host nests
	double Pa;						// probability of discovery of alien egg

	initializeParams(n, Pa);
	cuckooSearch(n, Pa);

	return 0;
}