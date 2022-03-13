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

void writeToFIles(vector<double>& X, vector<double>& Y, vector<double>& Z, vector<double>& Conv)
	{
	/*std::ifstream f1("fromY.csv");
	std::ifstream f2("fromIDFT.csv");
	std::ifstream f3("fromIFFT.csv");
	if (f1.good() || f2.good() || f3.good()) return;*/

	std::ofstream file1("x.txt");
	std::ofstream file2("y.txt");
	std::ofstream file3("z.txt");
	std::ofstream file4("conv.txt");

	doubleWrite(X, file1, "");
	doubleWrite(Y, file2, "");
	doubleWrite(Z, file3, "");
	doubleWrite(Conv, file4, "");
	
	file1.close();
	file2.close();
	file3.close();
	file4.close();
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
//x
void conv1(vector<double>& A, vector<double>& B, vector<double>& result)
	{
	//http://www.songho.ca/dsp/convolution/convolution.html
	int nconv = A.size() + B.size() - 1;
	int i, j, k, i1;

	for (i = B.size() - 1; i < A.size(); ++i)
		{
		result[i] = 0;                        
		for (j = i, k = 0; k < B.size(); --j, ++k)
			result[i] += A[j] * B[k];
		}

	for (i = 0; i < B.size() -1; ++i)
		{
		result[i] = 0;
		for (j = i, k = 0; j >= 0; --j, ++k)
			result[i] += A[j] * B[k];
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

	vector<double> XYConv(functionY.size() + functionZ.size());
	conv(functionY, functionZ, XYConv);
	writeToFIles(signalFunctionX, functionY, functionZ, XYConv);
#ifdef DBG
	for (int i=0; i< XYConv.size(); i++)
	std::cout << "XYConv " << XYConv[i] << " \n";
#endif
	}