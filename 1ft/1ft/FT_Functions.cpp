#pragma once
#include "FT_Functions.h"
#define DBG
#define F_FREQ
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
	double res = 0;
	size_t sizeAfterDFT = input.size();
		//https://scask.ru/c_book_r_cos.php?id=25
		for (size_t h = 0; h < sizeAfterDFT; h++)
			{
			auto ph = (2 * M_PI * h * numberOfSample_k) / sizeAfterDFT;
			res+= cos(ph) * input[h].real() - sin(ph) * input[h].imag();
			}
	res /= sizeAfterDFT;
	return res;
	}


// теряем знак
//double IDFT(vector<complex<double>>& input, int numberOfSample_k)
//	{
//	const auto complexI = std::complex<double>(0, 1);
//
//	size_t sizeAfterDFT = input.size();
//	std::complex<double> result;
//
//	for (size_t h = 0; h < sizeAfterDFT; h++)
//		{
//		//https://scask.ru/c_book_r_cos.php?id=25
//		result += exp((1. / sizeAfterDFT) * 2 * M_PI * h * numberOfSample_k * complexI) * input[h];
//		}
//	result /= sizeAfterDFT;
//	return std::abs(result);
//	}

// теряем знак
void FFT_AND_INVERS_FFT_TIME(vector<complex<double>>& a, bool invert) // по времени
	{
	int n = (int)a.size();
	if (n == 1)  return;

	vector<complex<double>> a0(n / 2), a1(n / 2);
	for (int i = 0, j = 0; i < n - 1; i += 2, ++j)
		{
		a0[j] = a[i];
		a1[j] = a[i + 1];
		}
	FFT_AND_INVERS_FFT_TIME(a0, invert);
	FFT_AND_INVERS_FFT_TIME(a1, invert);
	double ang = 2 * M_PI / n * (invert ? -1 : 1);
	complex<double> w(1), wn(cos(ang), sin(ang));
	for (int i = 0; i < n / 2; i++)
		{
		a[i] = a0[i] + w * a1[i];
		a[i + n / 2] = a0[i] - w * a1[i];
		if (invert)
			a[i] /= 2, a[i + n / 2] /= 2;
		w *= wn;
		}
	}

//херню считает
void FFT_AND_INVERS_FFT_FREQ(vector<complex<double>>& a, bool invert) // по времени
	{
	int n = (int)a.size();
	if (n == 1)  return;

	vector<complex<double>> b(n / 2);
	vector<complex<double>> c(n / 2);
	double ang = 2 * M_PI / n * (invert ? -1 : 1);
	complex<double> w(1), wn(cos(ang), sin(ang));

	for (int i = 0, j = 0; i < (n / 2); i += 1, ++j)
		{
		b[j] = a[i] + a[i + n / 2];
		c[j] = (a[i] - a[i + n / 2]) * w;
		w *= wn;
		} // k m p 
	FFT_AND_INVERS_FFT_FREQ(b, invert);
	FFT_AND_INVERS_FFT_FREQ(c, invert);
	for (int i = 0; i < n / 2; i++)
		{
		a[i] = b[i];
		a[i + n / 2] = c[i];
		}
	}
void fft_straight(vector<double> input, vector<complex<double>>& output)
	{
	for (int i = 0; i < input.size(); i++)
		{
		output[i] = (complex<double>(input[i]));
		}
#ifdef F_TIME
	cout << "Straight Time\n";
	FFT_AND_INVERS_FFT_TIME(output, 0);
#elif defined(F_FREQ)
	cout << "Straight Freq\n";
	FFT_AND_INVERS_FFT_FREQ(output, 0);
#endif 
#ifdef DBG
	for (int n=0; n<output.size(); n++)
	cout << "|" << n << "|" << output[n] << std::endl;
#endif 
	}

void fft_inverse(vector<complex<double>>& input, vector<double>& output)
	{
#if defined(F_TIME)
	cout << "Inverse Time\n";
	FFT_AND_INVERS_FFT_TIME(input, 1);
#elif defined(F_FREQ)
	cout << "Inverse Freq\n";
	FFT_AND_INVERS_FFT_FREQ(input, 1);
#endif 
	for (int i = 0; i < input.size(); i++) output[i] = std::abs(input[i]);
#ifdef DBG
	for (int i = 0; i < output.size(); i++)
	cout << "|" << i << "|" << output[i] << " \t>> module: " << abs(output[i]) << " phase: " << std::arg(output[i]) << std::endl;
#endif 
	}




void dft_straight(vector<double> input, vector<complex<double>>& output)
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