#pragma once
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
class PiNum
    {
        public:
        float _coeff;
        float num;
    PiNum(float coeff)
        {
        _coeff = coeff;
        num = PI * _coeff;
        }

    float PI = M_PI;

    float getPiNum()
        {
        return PI;
        }
    float getCoeff()
        {
        return _coeff;
        }
    float getNum()
        {
        return num;
        }

    };

