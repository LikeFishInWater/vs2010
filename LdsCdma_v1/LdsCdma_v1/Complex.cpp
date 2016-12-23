#include "Complex.h"


CComplex::CComplex(void)
	: Real(0)
	, Image(0)
{
}


CComplex::~CComplex(void)
{
}

CComplex::CComplex(double r, double i)
{
	Real = r;
	Image = i;
}

CComplex CComplex::operator+(CComplex & a)
{
	return CComplex(Real + a.Real, Image + a.Image);
}

CComplex CComplex::operator*(CComplex & a)
{
	CComplex temp;
	temp.Real = Real*a.Real - Image*a.Image;
	temp.Image = Image*a.Real + Real*a.Image;
	return temp;
}

CComplex CComplex::operator-(CComplex &a)
{
	return CComplex(Real - a.Real, Image - a.Image);
}
CComplex CComplex::operator/(CComplex &a)
{
	double temp = pow(a.Real, 2) + pow(a.Image, 2);
	return CComplex((Real*a.Real + Image*a.Image) / temp, (Image*a.Real - Real*a.Image) / temp);
}

double CComplex::Cabs()
{
	return sqrt(pow(Real, 2) + pow(Image, 2));
}
double CComplex::Cangle()
{
	double temp = atan(Image / Real);
	if (Image > 0)
		return temp;
	else
		return -temp;
}