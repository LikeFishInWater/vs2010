#pragma once
#include "math.h"
class CComplex
{
public:

	double Real;
	double Image;
public:
	CComplex(void);
	CComplex(double r, double i);
	~CComplex(void);
	CComplex operator+(CComplex & a);
	CComplex operator*(CComplex & a);
	
	CComplex operator-(CComplex &a);
	CComplex operator/(CComplex &a);
	double Cabs();
	double Cangle();
};

