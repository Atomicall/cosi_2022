//#define DBG
#include "main.h"
#include "FT_Functions.h"


#define N_O 64

using std::vector;
using std::complex;
using std::cout;
using namespace std::complex_literals;


float signal(float x)
	{
	return sin(sinCoeff * x) + cos(cosCoeff * x);
	}

bool compareDoubles(vector<double>& l, vector<double>& r)
	{
	if (l.size() != r.size()) return 0;
	if (l.empty() || r.empty()) return 0;
	auto it_r = r.begin();
	for (auto it_l = l.begin(); it_l != l.end(); it_l++, it_r++)
		{
		if (std::abs(*it_l - *it_r) > 0.1)
			{
			cout << "@@!= on " << *it_l << " - " << *it_r << std::endl;
			return 0;
			}
		
		}
	return 1;
	}


int main()
	{
	vector<double> signalFunctionX;
	vector<double> signalFunctionY;
	PiNum period = getPeriod();
	std::cout << "Func perios is " << period.getCoeff() << "Pi\n";
	float stepSize = period.getNum() / N_O;
	int index = 0;
	for (float i = 0; i < period.getNum(); i+=stepSize)
		{
		signalFunctionX.push_back(i);
		signalFunctionY.push_back( signal(i) );
#ifdef DBG
		std::cout << "For X " << signalFunctionX[index] << " Y: " << signalFunctionY[index] << " \n";
		index++;
#endif
		}
	vector<complex<double>> fromDFT;
	vector<double> fromInvertedDFT;

	vector<complex<double>> fromFFT(signalFunctionY.size());
	vector<double> fromInvertedFFT(signalFunctionY.size());
	dft_straight(signalFunctionY, fromDFT);
	dft_inverse(fromDFT, fromInvertedDFT);
	if (compareDoubles(signalFunctionY, fromInvertedDFT) == 1)
		{
		cout << "signalFunctionY == fromInvertedDFT\n";
		}
	fft_straight(signalFunctionY, fromFFT);
	fft_inverse(fromFFT, fromInvertedFFT);
	//if (compareDoubles(signalFunctionY, fromInvertedFFT) == 1) cout << "signalFunctionY == fromInvertedFFT\n";
	if (compareDoubles(signalFunctionY, fromInvertedFFT) == 1) cout << "signalFunctionY == fromInvertedFFT\n";
	writeToFIles(signalFunctionY, fromInvertedDFT, fromInvertedFFT, signalFunctionX);
	}


