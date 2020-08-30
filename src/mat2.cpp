#include "mat2.h"
#include <iostream>
#include <cmath>

using namespace std;

double& mat2::operator[](int t1)
{
	return x[t1];
}

mat2& mat2::operator=(const mat2& B)
{
	x[0] = B.x[0];
	x[1] = B.x[1];
	x[2] = B.x[2];
	x[3] = B.x[3];
	return *this;
}

mat2 mat2::operator+(const mat2& B)
{
	mat2 res;
	res.x[0] = x[0]+B.x[0];
	res.x[1] = x[1]+B.x[1];
	res.x[2] = x[2]+B.x[2];
	res.x[3] = x[3]+B.x[3];
	return res;
}

mat2 mat2::operator-(const mat2& B)
{
	mat2 res;
	res.x[0] = x[0]-B.x[0];
	res.x[1] = x[1]-B.x[1];
	res.x[2] = x[2]-B.x[2];
	res.x[3] = x[3]-B.x[3];
	return res;
}

mat2 mat2::operator*(const mat2& B)
{
	mat2 res;
	res.x[0] = x[0]*B.x[0]+x[1]*B.x[2];
	res.x[1] = x[0]*B.x[1]+x[1]*B.x[3];
	res.x[2] = x[2]*B.x[0]+x[3]*B.x[2];
	res.x[3] = x[2]*B.x[1]+x[3]*B.x[3];
	return res;
}

mat2 mat2::operator*(const double& B)
{
	mat2 res;
	res.x[0] = x[0]*B;
	res.x[1] = x[1]*B;
	res.x[2] = x[2]*B;
	res.x[3] = x[3]*B;
	return res;
}

mat2 mat2::operator/(const double& B)
{
	mat2 res;
	res.x[0] = x[0]/B;
	res.x[1] = x[1]/B;
	res.x[2] = x[2]/B;
	res.x[3] = x[3]/B;
	return res;
}


mat2 mat2::trans()
{
	mat2 res;
	res.x[0] = x[0];
	res.x[1] = x[2];
	res.x[2] = x[1];
	res.x[3] = x[3];
	return res;
}

void mat2::eig(mat2& vec, mat2& val)
{
	double e0 = (x[0]+x[3])/2;
	double delta = sqrt(pow((x[0]-x[3])/2,2)+x[1]*x[1]);
	val[0] = e0-delta;
	val[3] = e0+delta;
	vec[0] = x[3]-val[0];
	vec[2] = -x[1];
	vec[1] = -x[1];
	vec[3] = x[0]-val[3];
	double n1 = sqrt(vec[0]*vec[0]+vec[2]*vec[2]);
	double n2 = sqrt(vec[1]*vec[1]+vec[3]*vec[3]);
	vec[0] /= n1;
	vec[2] /= n1;
	vec[1] /= n2;
	vec[3] /= n2;
}


void mat2::print()
{
	cout<<x[0]<<'\t'<<x[1]<<'\n'<<x[2]<<'\t'<<x[3]<<endl;
}
