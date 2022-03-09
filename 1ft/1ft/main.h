#pragma once
#include "PiNum.h"
#include <fstream>
#include <vector>
#include <complex>
#include <iostream>

using std::vector;
using std::cout;
using std::complex;

#define sinCoeff 2
#define cosCoeff 7


unsigned int greatest_common_divisor(unsigned int a, unsigned int b)
	{
	if (a == b)
		return a;
	if (a > b)
		return greatest_common_divisor(a - b, b);
	return greatest_common_divisor(a, b - a);
	}
unsigned int least_common_multiple(unsigned int a, unsigned int b)
	{
	return (a * b) / greatest_common_divisor(a, b);
	}

PiNum getPeriod()
	{
	PiNum T1((float)2 / sinCoeff);
	PiNum T2((float)2 / cosCoeff);
	//return least_common_multiple(ceil(T1.getCoeff()), ceil(T2.getCoeff()));
	return PiNum(2);
	}

void writeToFIles(vector<double>& fromY, vector<double>& fromIDFT, vector<double>& fromIFFT, vector<double>& fromX)
	{
	std::ifstream f1("fromY.csv");
	std::ifstream f2("fromIDFT.csv");
	std::ifstream f3("fromIFFT.csv");
	if (f1.good() || f2.good() || f3.good()) return;

	std::ofstream file1("fromY.csv");
	std::ofstream file2("fromIDFT.csv");
	std::ofstream file3("fromIFFT.csv");
	file1 << "fromX;fromY\n";
	for (size_t i = 0; i < fromY.size(); i++)
		{
		file1 << std::fixed << fromX[i] << ";" << fromY[i] << "\n";
		}

	file2 << "fromX;fromIDFT\n";
	for (size_t i = 0; i < fromIDFT.size(); i++)
		{
		file2 << std::fixed << fromX[i] << ";" << fromIDFT[i] << "\n";
		}

	file3 << "fromX;fromIFFT\n";
	for (size_t i = 0; i < fromIFFT.size(); i++)
		{
		file3 << std::fixed << fromX[i] << ";" << fromIFFT[i] << "\n";
		}
	file1.close();
	file2.close();
	file3.close();
	}


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

void writeAmplPhaseToFIles(vector<complex<double>>& dft, vector<complex<double>>& fft)
	{
	std::ifstream f1("amplDFT.csv");
	std::ifstream f2("phaseDFT.csv");
	std::ifstream f3("amplFFT.csv");
	std::ifstream f4("phaseFFT.csv");
	if (f1.good() || f2.good() || f3.good() || f4.good()) return;

	std::ofstream file1("amplDFT.csv");
	std::ofstream file2("phaseDFT.csv");
	std::ofstream file3("amplFFT.csv");
	std::ofstream file4("phaseFFT.csv");

	writeComplexVector(dft, file1, file2, "DFT");
	writeComplexVector(fft, file3, file4, "FFT");
	file1.close();
	file2.close();
	file3.close();
	file4.close();
	}


