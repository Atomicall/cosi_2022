#pragma once
#include "PiNum.h"

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
