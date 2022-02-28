#include "main.h"
#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>

#define N 64
#define sinCoeff 2
#define cosCoeff 7



void dft_straight(float* samples, int& samples_amount, float* result);
void dft_reverse(float* samples, int& samples_amount, float* result);
void fft_straight(float* samples, int& samples_amount, float* result);
void fft_reverse(float* samples, int& samples_amount, float* result);
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

PiNum getPeriod()
	{
	PiNum T1((float)2 / sinCoeff);
	PiNum T2((float)2/cosCoeff);
	//return least_common_multiple(ceil(T1.getCoeff()), ceil(T2.getCoeff()));
	return PiNum(2);
	}

int main()
	{
	float* signalFunctionX = new float[N];
	float* signalFunctionY = new float[N];
	PiNum period = getPeriod();
	std::cout << "Func perios is " << period.getCoeff() << "Pi\n";
	float stepSize = period.getNum() / N;
	int index = 0;
	for (float i = 0; i < period.getNum(); i+=stepSize)
		{
		signalFunctionX[index] = i;
		signalFunctionY[index] = signal(i);
		//std::cout << "For X " << signalFunctionX[index] << " Y: " << signalFunctionY[index] << " \n";
		index++;
		}



	}

void dft_straight(float* samples, int& samples_amount, float* result)
	{
	}

void dft_reverse(float* samples, int& samples_amount, float* result)
	{
	}

void fft_straight(float* samples, int& samples_amount, float* result)
	{
	}

void fft_reverse(float* samples, int& samples_amount, float* result)
	{
	}
