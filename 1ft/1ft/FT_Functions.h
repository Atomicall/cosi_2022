#pragma once
#include <vector>
#include <complex>
#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>


using std::vector;
using std::complex;
using std::cout;
using namespace std::complex_literals;

void dft_straight(vector<double> input, vector<complex<double>>& output);
void dft_inverse(vector<complex<double>>& input, vector<double>& output);
void fft_straight(vector<double> input, vector<complex<double>>& output);
void fft_inverse(vector<complex<double>>& input, vector<double>& output);

std::tuple<complex<double>, std::pair<int, int>> DFT(int numberOfSample_k, vector<double>& samples);
std::tuple<double, std::pair<int, int>> IDFT(vector<complex<double>>& input, int numberOfSample_k);
void FFT_AND_INVERS_FFT_TIME(vector<complex<double>>& a, bool invert);
void FFT_AND_INVERS_FFT_FREQ(vector<complex<double>>& a, bool invert);
