// ati AGE_Tema0.cpp : Defines the entry point for the console application.
#include <stdafx.h>
#include <iostream>
#include <ctime> 
#include <cstdlib>
#include <cmath>
#include <math.h>
#define sizeOfSample 10
#define numberOfRuns 100000
#define maxNoOfParam 3
#define MIN_VAL 9999999.0
#define MAX_VAL 0.0
#define PI 3.1415926535897
#define DMAX 1000

using namespace std;

double TEMP = 100;
double t_MINUS = 0.98;
int numberOfBits;
int noOfParameters;
enum Name { DeJong, Schwefels, Rastrigins } name;

struct Interval {
	double lowValue;
	double highValue;
}param1, param2, param3, values[maxNoOfParam];

int Calculate_noOfBits(double a, double b) {

	//10^d, d = 2;
	double N = (b - a) * 100;
	double n = log2(N);
	int result = n;
	if (n == result)
		return result;
	return result + 1;
}

double randomValue(double lowValue, double highValue) {
	return lowValue + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (highValue - lowValue)));
}

double DeJong1(double sample[sizeOfSample]) {
	double sum = 0;
	for (int index = 0; index < sizeOfSample; index++) {
		sum += sample[index] * sample[index];
	}
	return sum;
}

double Schwefel(double sample[sizeOfSample]) {
	double sum = 0;
	for (int index = 0; index < sizeOfSample; index++) {
		sum += (-1) * sample[index] * sin(sqrt(abs(sample[index])));
	}
	return sum;
}

double Rastrigin(double sample[sizeOfSample]) {
	double sum = 0;
	for (int index = 0; index < sizeOfSample; index++) {
		sum += (sample[index] * sample[index] - 10 * cos(2 * PI*sample[index]));
	}
	return sum;
}

double CalculateResult(Name functionName, double sample[sizeOfSample]) {

	double result;
	switch (functionName) {
	case DeJong:
		result = DeJong1(sample);
		break;
	case Schwefels:
		result = Schwefel(sample);
		break;
	case Rastrigins:
		result = 10 * sizeOfSample + Rastrigin(sample);
		//result = Rastrigin(sample);
		break;
	default:
		break;
	}
	return result;
}

double ToBase10(Name name, int vector[], int length) {

	double number = 0;
	for (int index = length - 1; index >= 1; index--) {
		if (vector[index] == 1) {
			number += pow(2, length - index - 1);
		}
	}
	return number / 100;
}

double Random_Neighbour(Name name, int vector[], int length, double sample[sizeOfSample], int random, int rand) {

	double min = MIN_VAL, result = 0;
	int position = 0;

	vector[rand] = 1 - vector[rand];
	sample[random] = ToBase10(name, vector, numberOfBits);
	if (vector[0] == 1)
		sample[random] *= -1;
	result = CalculateResult(name, sample);

	return result;
}

void SimulatedAnnealing(Name functionName, int noOfParam, Interval values[]) {

	//int M[sizeOfSample][numberOfBits];
	int runs, i = 0;
	int rows = sizeOfSample, cols = numberOfBits;
	int **M;
	M = new int*[rows];
	for (int i = 0; i<rows; i++) {
		M[i] = new int[cols];
	}
	double sample[sizeOfSample], minSample[sizeOfSample];
	double minimum, localMin, globalMin = MIN_VAL, result = 0;
	int rValue;
	bool ok;
	for (int parameters = 0; parameters < noOfParam; parameters++) {

		minimum = MIN_VAL; globalMin = MIN_VAL;

		for (int index = 0; index < sizeOfSample; index++) {
			for (int bit = 0; bit < numberOfBits; bit++) {
				M[index][bit] = rand() % 2;
			}
		}//am generat o matrice random de biti

		for (int index = 0; index < sizeOfSample; index++) {
			sample[index] = ToBase10(functionName, M[index], numberOfBits);
			if (M[index][0] == 1)
				sample[index] *= -1;
		}
		minimum = CalculateResult(functionName, sample); //aflam valoarea Vc

		while ( TEMP > 0.00001){

			for (runs = 0; runs < numberOfRuns; runs++) {
				
				ok = false;
				rValue = rand() % sizeOfSample; //selectez un vecin random
				int rPos = rand() % numberOfBits + 1;

				localMin = Random_Neighbour(functionName, M[rValue], numberOfBits, sample, rValue, rPos);
				if (localMin < minimum) {
					minimum = localMin;
					for (int index = 0; index < sizeOfSample; index++)
						minSample[index] = sample[index];
					M[rValue][rPos] = 1 - M[rValue][rPos];
				}
				else {
					if (randomValue(0, 0.999999999) < exp(-abs(localMin - minimum) / TEMP)) {
						minimum = localMin;
						for (int index = 0; index < sizeOfSample; index++)
							minSample[index] = sample[index];
						M[rValue][rPos] = 1 - M[rValue][rPos];
					}
				}
			}

			if (minimum < globalMin) {
				globalMin = minimum;
				cout << globalMin << '\n';
			}
			TEMP *= t_MINUS;

		}
		for (int parameters = 0; parameters < noOfParam; parameters++) {
			cout << "X_" << parameters + 1 << ": \n";
			for (int index = 0; index < sizeOfSample; index++) {
				cout << minSample[index] << " ";
			}
		}
		cout << "\nY_MIN: " << globalMin << '\n';
		cout << "TEMP: " << TEMP << '\n';
		cout << runs;
	}
}

void DeJong_Init() {
	name = DeJong;
	param1.lowValue = -5.12;
	param1.highValue = 5.12;
	values[0] = param1;
	noOfParameters = 1;

	numberOfBits = Calculate_noOfBits(param1.lowValue, param1.highValue);
}

void Schwefels_Init() {
	name = Schwefels;
	param1.lowValue = -500.0;
	param1.highValue = 500.0;
	values[0] = param1;
	noOfParameters = 1;
	numberOfBits = Calculate_noOfBits(param1.lowValue, param1.highValue);
}

void Rastrigins_Init() {
	name = Rastrigins;
	param1.lowValue = -5.12;
	param1.highValue = 5.12;
	values[0] = param1;
	noOfParameters = 1;
	numberOfBits = Calculate_noOfBits(param1.lowValue, param1.highValue);
}

int main() {

	srand((unsigned int)time(NULL));
	//DeJong_Init();
	Schwefels_Init();
	//Rastrigins_Init();
	SimulatedAnnealing(name, noOfParameters, values);
	return 0;
}
