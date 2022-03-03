#include "main.h"
#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <vector>
#include <complex>

//#define DBG

#define N_O 64
#define sinCoeff 2
#define cosCoeff 7


using std::vector;
using std::complex;
using std::cout;
using namespace std::complex_literals;


void dft_straight(vector<double>& input, vector<complex<double>>& output);
void dft_inverse(vector<complex<double>>& input, vector<double>& output);
void fft_straight(vector<double>& input, vector<complex<double>>& output);
void fft_inverse(vector<complex<double>>& input, vector<double>& output);


float signal(float x)
	{
	return sin(sinCoeff * x) + cos(cosCoeff * x);
	}
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
bool compareDoubles(vector<double>& l, vector<double>& r)
	{
	if (l.size() != r.size()) return 0;
	auto it_r = r.begin();
	for (auto it_l = l.begin(); it_l != l.end(); it_l++, it_r++)
		{
		if (*it_l - *it_r > 0.1)
			{
			cout << "!= on " << *it_l << " - " << *it_r << std::endl;
			return 0;
			}
		}
	return 1;
	}
PiNum getPeriod()
	{
	PiNum T1((float)2 / sinCoeff);
	PiNum T2((float)2/cosCoeff);
	//return least_common_multiple(ceil(T1.getCoeff()), ceil(T2.getCoeff()));
	return PiNum(2);
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
		//std::cout << "For X " << signalFunctionX[index] << " Y: " << signalFunctionY[index] << " \n";
		index++;
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
	if (compareDoubles(signalFunctionY, fromInvertedFFT)) cout << "signalFunctionY == fromInvertedFFT\n";


	}






complex<double> DFT(int numberOfSample_k, vector<double>& samples)
	{
	//https://ru.dsplib.org/dspl/group___d_f_t___g_r_o_u_p.html

	double a = 0;
	double b = 0;
	auto samplesAmount = samples.size(); //64

	complex<double> teeemp(0, 0);
	for (int n = 0; n < samplesAmount; n++)
		{
		// - i * (2 PI numberOfSample_k cycleI )/(OverallNumberOfSamlesInInput )
		//numberOfSample_k - номер отсчета из ВХОДНОЙ послед для котороко ищем C
		// 
		a += cos((2 * M_PI * numberOfSample_k * n) / samplesAmount) * samples[n]; // Действительная
		b += -sin((2 * M_PI * numberOfSample_k * n) / samplesAmount) * samples[n]; // Мнимая

		/*complex<double> aa(0, ((2 * M_PI * numberOfSample_k * n) / (samplesAmount)));
		teeemp += samples[n] * aa;*/

		}
	complex<double> Ck(a, b);
	return Ck;

	//return teeemp;
	}
double IDFT(vector<complex<double>>& input, int numberOfSample_k)
	{

	const auto complexI = std::complex<double>(0, 1);

	size_t sizeAfterDFT = input.size();
	std::complex<double> result;

	for (size_t h = 0; h < sizeAfterDFT; h++)
		{
		//https://scask.ru/c_book_r_cos.php?id=25
		result += exp((1. / sizeAfterDFT) * 2 * M_PI * h * numberOfSample_k * complexI) * input[h];
		}
	result /= sizeAfterDFT;
	return std::abs(result);
	}

void dft_straight(vector<double>& input, vector<complex<double>>& output)
	{
	for (int i = 0; i < input.size(); i++)
		{
		complex<double> resultForIElement = DFT(i, input); // k-й результат до k-го элемента из input // Ck-е
		//https://scask.ru/a_book_tec.php?id=149
		output.push_back(resultForIElement);
		
#ifdef DBG
		cout << "|" << i << "|" << output[i] << " \t>> module: " << abs(output[i]) << " phase: " << arg(output[i]) << std::endl;
#endif 
		}
	}

void dft_inverse(vector<complex<double>>& input, vector<double>& output)
	{
	//https://ru.dsplib.org/dspl/group___d_f_t___g_r_o_u_p.html 
	for (int n = 0; n < input.size(); n++)
		{
		output.push_back(IDFT(input, n));
#ifdef DBG
		cout << "|" << n << "|" << output[n] << std::endl;
#endif 
		}
	}


void FFT_AND_INVERS_FFT_TIME(vector<complex<double>>& a, bool invert) // по времени
	{
	int n = (int)a.size();
	if (n == 1)  return;

	vector<complex<double>> a0(n / 2), a1(n / 2);
	for (int i = 0, j = 0; i < n-1; i += 2, ++j)
		{
		a0[j] = a[i];
		a1[j] = a[i + 1];
		}
	FFT_AND_INVERS_FFT_TIME(a0, invert);
	FFT_AND_INVERS_FFT_TIME(a1, invert);
	double ang = 2 * M_PI / n * (invert ? -1 : 1);
	complex<double> w(1), wn(cos(ang), sin(ang));
	for (int i = 0; i < n / 2; ++i)
		{
		a[i] = a0[i] + w * a1[i];
		a[i + n / 2] = a0[i] - w * a1[i];
		if (invert)
			a[i] /= 2, a[i + n / 2] /= 2;
		w *= wn;
		}
	}

void fft_straight(vector<double>& input, vector<complex<double>>& output)
	{
	for (int i = 0; i < input.size(); i++)
		{
		output[i] = (complex<double>(input[i]));
		}
	FFT_AND_INVERS_FFT_TIME(output, 0);


	}

void fft_inverse(vector<complex<double>>& input, vector<double>& output)
	{
	FFT_AND_INVERS_FFT_TIME(input, 1);
	for (int i = 0; i < input.size(); i++) output[i] = std::abs(input[i]);
	}


