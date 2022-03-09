#pragma once
#include "FT_Functions.h"
#include <tuple>
#define DBG
#define F_FREQ
std::tuple<complex<double>, std::pair<int, int>> DFT(int numberOfSample_k, vector<double>& samples)
	{
	//https://ru.dsplib.org/dspl/group___d_f_t___g_r_o_u_p.html

	static int adds = 0, muls = 0, maxsize = 0;
	

	double a = 0;
	double b = 0;
	auto samplesAmount = samples.size(); //64
	if (maxsize <= samplesAmount)
		{
		maxsize = samplesAmount;
		adds = muls = 0;
		}

	complex<double> teeemp(0, 0);
	for (int n = 0; n < samplesAmount; n++)
		{
		// - i * (2 PI numberOfSample_k cycleI )/(OverallNumberOfSamlesInInput )
		//numberOfSample_k - номер отсчета из ВХОДНОЙ послед для котороко ищем C
		// 
		a += cos((2 * M_PI * numberOfSample_k * n) / samplesAmount) * samples[n]; // Действительная
		b += -sin((2 * M_PI * numberOfSample_k * n) / samplesAmount) * samples[n]; // Мнимая
		adds += 1;
		muls += 5;
		/*complex<double> aa(0, ((2 * M_PI * numberOfSample_k * n) / (samplesAmount)));
		teeemp += samples[n] * aa;*/

		}
	complex<double> Ck(a, b);
	auto p = std::make_pair<int, int>((int)adds, (int)muls);
	return std::make_tuple(Ck, p);
//#ifdef DBG
//	std::cout << "DFT: adds:" << adds << " muls" << muls << " \n";
//#endif 
	//return Ck;
	//return teeemp;
	}
std::tuple<double, std::pair<int, int>> IDFT(vector<complex<double>>& input, int numberOfSample_k)
	{
	static int adds = 0, muls = 0, maxsize = 0;
	double res = 0;
	size_t sizeAfterDFT = input.size();
	if (maxsize <= sizeAfterDFT)
		{
		maxsize = sizeAfterDFT;
		adds = muls = 0;
		}
		//https://scask.ru/c_book_r_cos.php?id=25
		for (size_t h = 0; h < sizeAfterDFT; h++)
			{
			auto ph = (2 * M_PI * h * numberOfSample_k) / sizeAfterDFT;
			res+= cos(ph) * input[h].real() - sin(ph) * input[h].imag();
			adds += 1;
			muls += 5;
			}
	res /= sizeAfterDFT;
	muls++;
	auto p = std::make_pair<int, int>((int)adds, (int)muls);
	return std::make_tuple(res, p);
//#ifdef DBG
//	std::cout << "IDFT: adds:" << adds << " muls" << muls << " \n";
//#endif 

	//return res;
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

////херню считает
//void FFT_AND_INVERS_FFT_FREQ1(vector<complex<double>>& a, bool invert) // по времени
//	{
//	int n = (int)a.size();
//	if (n == 1)  return;
//
//	vector<complex<double>> b(n / 2);
//	vector<complex<double>> c(n / 2);
//	double ang = 2 * M_PI / n * (invert ? -1 : 1);
//	complex<double> w(1), wn(cos(ang), sin(ang));
//
//	for (int i = 0, j = 0; i < (n / 2); i += 1, ++j)
//		{
//		b[j] = a[i] + a[i + n / 2];
//		c[j] = (a[i] - a[i + n / 2]) * w;
//		w *= wn;
//		} // k m p 
//	FFT_AND_INVERS_FFT_FREQ1(b, invert);
//	FFT_AND_INVERS_FFT_FREQ1(c, invert);
//	for (int i = 0; i < n / 2; i++)
//		{
//		a[i] = b[i];
//		a[i + n / 2] = c[i];
//		}
//	}


void FFT_AND_INVERS_FFT_FREQ(vector<complex<double>>& a, bool invert) // по частоте
	{
	static int adds = 0, muls = 0, maxsize = 0;
	int n = (int)a.size();
	if (maxsize <= n)
		{
		maxsize = n;
		adds = muls = 0;
		}
	if (n == 1)  return;

	vector<complex<double>> y(n / 2);
	vector<complex<double>> z(n / 2);

	for (size_t i = 0; i < n/2; i++)
		{
		adds+=2;
		muls+=2;
		y[i] = a[i] + a[i + (n / 2)];
		z[i] = (a[i] - a[i + (n / 2)]) * std::polar(1.0, (invert ? -1 : 1) * -2 * M_PI * i / (double)n);
		}
	FFT_AND_INVERS_FFT_FREQ(y, invert);
	FFT_AND_INVERS_FFT_FREQ(z, invert);

	for (int i = 0; i < n / 2; i++)
		{
		a[2* i] = y[i];
		a[2*i + 1] = z[i];
		if (invert)
			a[2 * i] /= 2, a[2 * i + 1] /= 2, muls += 2;
		}
#ifdef DBG
	std::cout << "FFT_AND_INVERS_FFT_FREQ: adds:" << adds << " muls" << muls << " on invert = " << std::to_string(invert) << " \n";
#endif 
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
	const auto complexI = std::complex<double>(0, 1);
	for (int i = 0; i < input.size(); i++) { 
		output[i] = input[i].real();
		};

#ifdef DBG
	for (int i = 0; i < output.size(); i++)
	cout << "|" << i << "|" << output[i] << " \t>> module: " << abs(output[i]) << " phase: " << std::arg(output[i]) << std::endl;
#endif 
	}

void dft_straight(vector<double> input, vector<complex<double>>& output)
	{
	int adds = 0, muls = 0;
	for (int i = 0; i < input.size(); i++)
		{
		complex<double> resultForIElement;
		auto res = DFT(i, input); // k-й результат до k-го элемента из input // Ck-е
		resultForIElement = std::get<0>(res);
		adds += std::get<1>(res).first;
		muls += std::get<1>(res).second;
		//https://scask.ru/a_book_tec.php?id=149
		output.push_back(resultForIElement);

#ifdef DBG
		cout << "|" << i << "|" << output[i] << " \t>> module: " << abs(output[i]) << " phase: " << arg(output[i]) << std::endl;
#endif 
		}
#ifdef DBG
	std::cout << "dft_straight: adds:" << adds << " muls " << muls << " \n";
#endif
	}

void dft_inverse(vector<complex<double>>& input, vector<double>& output)
	{
	int adds = 0, muls = 0;
	//https://ru.dsplib.org/dspl/group___d_f_t___g_r_o_u_p.html 
	for (int n = 0; n < input.size(); n++)
		{
		auto result = IDFT(input, n);
		output.push_back(std::get<0>(result));
		adds += std::get<1>(result).first;
		muls += std::get<1>(result).second;
#ifdef DBG
		cout << "|" << n << "|" << output[n] << std::endl;
#endif 
		}
#ifdef DBG
	std::cout << "dft_inverse: adds:" << adds << " muls " << muls << " \n";
#endif
	}