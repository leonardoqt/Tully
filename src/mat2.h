#ifndef __MAT2__
#define __MAT2__

class mat2
{
public:
	double x[4];
	
	double& operator[](int);
	mat2& operator=(const mat2&);
	mat2 operator+(const mat2&);
	mat2 operator-(const mat2&);
	mat2 operator*(const mat2&);
	mat2 operator*(const double&);
	mat2 operator/(const double&);
	
	mat2 trans();
	void eig(mat2& vec, mat2& val);

	void print();
};
#endif
