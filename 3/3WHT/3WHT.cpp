
#define DBG
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <iostream>
//y=sin(2x)+cos(7x)

#define N_O 63

using std::vector;
using std::complex;
using std::cout;
using namespace std::complex_literals;
#define sinCoeff 2
#define cosCoeff 7
float signal(float x)
    {
    return sin(sinCoeff * x) + cos(cosCoeff * x);
    }

void fwht(vector<complex<double>>& input) //input <=output
	{
	int h = 1;
	while (h < input.size())
		{
		for (int i = 0; i < input.size(); i += h * 2)
			{
			for (int j = i; j < i + h; j++)
				{
				auto x = input[j];
				auto y = input[j + h];
				input[j] = x + y;
				input[j + h] = x - y;
				}
			}
		h *= 2;
		}
	}

//def fwht(a)->None:
//"""In-place Fast Walsh–Hadamard Transform of array a."""
//h = 1
//while h < len(a) :
//    for i in range(0, len(a), h * 2) :
//        for j in range(i, i + h) :
//            x = a[j]
//            y = a[j + h]
//            a[j] = x + y
//            a[j + h] = x - y
//            h *= 2

void writeComplexVector(vector<complex<double>>& what, std::ofstream& whereAmpl, std::ofstream& wherePhase, std::string header)
	{
	whereAmpl << "N;" << header << "\n";
	wherePhase << "N;" << header << "\n";
	for (size_t i = 0; i < what.size(); i++)
		{
		whereAmpl << i << ";" << std::abs(what[i]) << "\n";
		wherePhase << i << ";" << std::arg(what[i]) << "\n";
		}
	}

void writeAmplPhaseToFIles(vector<complex<double>>& fwdht)
	{
	std::ifstream f1("FWHT_A.csv");
	std::ifstream f2("FWHT_P.csv");
	if (f1.good() || f2.good()) return;
	std::ofstream file1("FWHT_A.csv");
	std::ofstream file2("FWHT_P.csv");
	writeComplexVector(fwdht, file1, file2, "FWHT");
	file1.close();
	file2.close();
	}

int main()
	{
	vector<double> signalFunctionX;
	vector<double> signalFunctionY;
	float stepSize = (2 * M_PI )/ N_O;
	int index = 0;
	for (float i = 0; i < 2 * M_PI; i += stepSize)
		{
		signalFunctionX.push_back(i);
		signalFunctionY.push_back(signal(i));
		}
#ifdef DBG
		for (int i = 0; i < signalFunctionX.size(); i++)
		std::cout << "For X " << signalFunctionX[i] << " Y: " << signalFunctionY[i] << " \n";
#endif
		vector<complex<double>> WHTResult(signalFunctionY.size());
		for (size_t i = 0; i < WHTResult.size(); i++)
			{
			WHTResult[i] = complex<double>(signalFunctionY[i]);
			}
		fwht(WHTResult);
#ifdef DBG
		for (int i=0; i< WHTResult.size(); i++)
		std::cout << "WHT mod" << std::abs(WHTResult[i]) << "WHT phase" << std::arg(WHTResult[i]) << " \n";
#endif
		writeAmplPhaseToFIles(WHTResult);

	}

