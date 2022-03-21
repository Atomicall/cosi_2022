#include <fstream>
#include <vector>
#include <complex>
#include <iostream>
#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include "misc.h"
#define DBG
//4 вариант
// y=sin(2x)
// z =cos(7x)
// N 16

using std::vector;
using std::cout;
using std::complex;

#define N_O 16
#define sinCoeff 2.0
#define cosCoeff 7.0

float signalY(float x)
	{
	return sin(sinCoeff * x);
	}

float signalZ(float x)
	{
	return cos(cosCoeff * x);
	}

void doubleWrite(vector<double>& what, std::ofstream& where, std::string head)
	{
	where << head << "\n";
	for (size_t i = 0; i < what.size(); i++)
		{
		where << std::fixed << what[i] << " ";
		}
	}
void writeToFIles(vector<double>& X, vector<double>& Y, vector<double>& Z, vector<double>& Conv, vector<double>& Corr)
	{
	/*std::ifstream f1("fromY.csv");
	std::ifstream f2("fromIDFT.csv");
	std::ifstream f3("fromIFFT.csv");
	if (f1.good() || f2.good() || f3.good()) return;*/

	std::ofstream file1("x.txt");
	std::ofstream file2("y.txt");
	std::ofstream file3("z.txt");
	std::ofstream file4("conv.txt");
	std::ofstream file5("corr.txt");

	doubleWrite(X, file1, "");
	doubleWrite(Y, file2, "");
	doubleWrite(Z, file3, "");
	doubleWrite(Conv, file4, "");
	doubleWrite(Corr, file5, "");
	
	file1.close();
	file2.close();
	file3.close();
	file4.close();
	file5.close();
	}

void conv(vector<double>& A, vector<double>& B, vector<double>& result)
	{
	//https://toto-share.com/2011/11/cc-convolution-source-code/
	int nconv = A.size() + B.size() - 1;
	int i, j, i1;
	float tmp;
	for (i = 0; i < nconv; i++)
		{
		i1 = i;
		tmp = 0.0;
		for (j = 0; j < B.size(); j++)
			{
			if (i1 >= 0 && i1 < A.size())
				{
				tmp = tmp + (A[i1] * B[j]);
				}
			i1 = i1 - 1;
			result[i] = tmp;
			}
		}
	}

void corr(vector<double>& A, vector<double>& B, vector<double>& result)
	{
	vector<double> a;
	vector<double> b;
	vector<double> bTemp(B);
	size_t N = A.size();
	size_t M = B.size();
	a.resize(A.size() + M - 1, 0);
	b.resize(B.size() + N - 1, 0);
	for (size_t i = M - 1, j = 0; i < a.size(); i++, j++)
		{
		a[i] = A[j];
		}
	std::reverse(bTemp.begin(), bTemp.end());
	for (size_t i = 0; i < bTemp.size(); i++)
		{
		b[i] = bTemp[i];
		}
	double res = 0.0;
	for (size_t m = 0; m < N + M - 1; m++)
		{
		for (size_t i = 0; i < N + M - 1; i++)
			{
			res += a[i] * b[i];
			}
		result[m] = res;
		res = 0.0;
		for (size_t j = N+M-2; j > m; j--)
			{
			b[j] = b[j - 1];
			}
		b[m] = 0;
		}
	}


void fft_time(vector<complex<double>>& arr)
	{
	const size_t N = arr.size();
	if (N <= 1) return;

	vector<complex<double>> even = arr[std::slice(0, N / 2, 2)];
	vector<complex<double>> odd = arr[std::slice(1, N / 2, 2)];

	fft_time(even);
	fft_time(odd);

	for (int j = 0; j < N / 2; j++)
		{
		complex<double> w = std::polar(1.0, -2 * M_PI * j / N) * odd[j];
		arr[j] = even[j] + w;
		arr[j + N / 2] = even[j] - w;
		}
	}

void ifft(vector<complex<double>>& x)
	{
	// conjugate the complex numbers
	for (auto it = x.begin(); it != x.end() ; it++)
		{
		*it = std::conj((*it));
		}
	// forward fft
	fft_time(x);
	// conjugate the complex numbers again
	for (auto it = x.begin(); it != x.end(); it++)
		{
		*it = std::conj((*it));
		}
	// scale the numbers
	x /= x.size();
	for (auto it = x.begin(); it != x.end(); it++)
		{
		*it = *it / x.size();
		}
	}


int main()
	{
	vector<double> signalFunctionX;
	vector<double> functionY;
	vector<double> functionZ;
	float stepSize = (2 * M_PI )/ N_O;
	int index = 0;
	for (float i = 0; i < 2* M_PI; i += stepSize)
		{
		signalFunctionX.push_back(i);
		functionY.push_back(signalY(i));
		functionZ.push_back(signalZ(i));
#ifdef DBG
		std::cout << "For X " << signalFunctionX[index] << " Y: " << functionY[index] << " Z: " << functionZ[index] << " \n";
		index++;
#endif
		}

	vector<double> XYConv(functionY.size() + functionZ.size()-1);
	vector<double> XYCorr(functionY.size() + functionZ.size()-1);

	conv(functionY, functionZ, XYConv);
	corr(functionY, functionZ, XYCorr);

	writeToFIles(signalFunctionX, functionY, functionZ, XYConv, XYCorr);
#ifdef DBG
	std::cout << "XYConv " << "\n";
	for (int i = 0; i < XYConv.size(); i++)
		{
		std::cout << XYConv[i] << " ; ";
		if (!(i%3))	std::cout << "\n";
		}
	std::cout << "\n";
#endif

#ifdef DBG
	std::cout << "XYCorrelatoin " << "\n";
	for (int i = 0; i < XYCorr.size(); i++)
		{
		std::cout << XYCorr[i] << " ; ";
		if (!(i % 3))	std::cout << "\n";
		}
	std::cout << "\n";
#endif
	}