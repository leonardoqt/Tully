#ifndef __POT__
#define __POT__

#include "mat2.h"
#include <cmath>

mat2 H1(double x);
mat2 H2(double x);
mat2 H3(double x);


mat2 H1(double x)
{
	mat2 HH;
	double A=0.01,B=1.6,C=0.005,D=1.0;
	if (x >= 0)
		HH[0] = A*(1-std::exp(-B*x));
	else
		HH[0] = -A*(1-std::exp(B*x));
	HH[2] = HH[1] = C*std::exp(-D*x*x);
	HH[3] = -HH[0];
	return HH;
}

mat2 H2(double x)
{
	mat2 HH;
	double A=0.1,B=0.28,E0=0.05,C=0.015,D=0.06;
	HH[0] = 0;
	HH[2] = HH[1] = C*std::exp(-D*x*x);
	HH[3] = -A*std::exp(-B*x*x) + E0;
	return HH;
}

mat2 H3(double x)
{
	mat2 HH;
	double A = 6.0e-4,B=0.1,C=0.9;
	HH[0] = A;
	HH[3] = -A;
	if (x <=0)
		HH[2] = HH[1] = B*std::exp(C*x);
	else
		HH[2] = HH[1] = B*(2-std::exp(-C*x));
	return HH;
}

#endif
