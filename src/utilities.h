#ifndef __UTIT__
#define __UTIT__

#include "mat2.h"
#include "pot.h"
#include <cmath>

mat2 grad_E(double x, mat2 (*Hp)(double));
mat2 DD(double x, mat2 (*Hp)(double));

void x_v_rk4(double dT, double& x, double& v, double m, int s, mat2 (*Hp)(double));
void rho_rk4(double dT, mat2& rho, double& rho12im, mat2& rho_dot, double v, mat2 DD, mat2 EE);
//==========================================

mat2 grad_E(double x, mat2 (*Hp)(double))
{
	double dx = 1.0e-6;
	mat2 res,vec,val1,val2,HH;
	HH = Hp(x+dx);
	HH.eig(vec,val2);
	HH = Hp(x-dx);
	HH.eig(vec,val1);
	res = (val2-val1)/(2*dx);
	return res;
}

mat2 DD(double x, mat2 (*Hp)(double))
{
	double dx = 1.0e-6;
	mat2 HH, res,vec,val;
	HH = Hp(x);
	HH.eig(vec,val);
	HH = (Hp(x+dx) - Hp(x-dx))/(2*dx);
	res = (vec.trans()*HH)*vec;
	res[0] = res[3] = 0;
	res[1] /= val[3]-val[0];
	res[2] /= val[0]-val[3];
	return res;
}

void x_v_rk4(double dT, double& x, double& v, double m, int s, mat2 (*Hp)(double))
{
	double xk1,xk2,xk3,xk4;
	double vk1,vk2,vk3,vk4;
	mat2 grad;
	// k1
	xk1 = v;
	grad = grad_E(x,Hp);
	vk1 = -grad[s]/m;
	// k2
	xk2 = v + vk1*dT/2;
	grad = grad_E(x+xk1*dT/2,Hp);
	vk2 = -grad[s]/m;
	// k3
	xk3 = v + vk2*dT/2;
	grad = grad_E(x+xk2*dT/2,Hp);
	vk3 = -grad[s]/m;
	// k4
	xk4 = v + vk3*dT;
	grad = grad_E(x+xk3*dT,Hp);
	vk4 = -grad[s]/m;
	// final
	x = x + (xk1+2*xk2+2*xk3+xk4)*dT/6;
	v = v + (vk1+2*vk2+2*vk3+vk4)*dT/6;
}

void rho_rk4(double dT, mat2& rho, double& rho12im, mat2& rho_dot, double v, mat2 DD, mat2 EE)
{
	double r11k1,r11k2,r11k3,r11k4;
	double r12k1,r12k2,r12k3,r12k4;
	double i12k1,i12k2,i12k3,i12k4;

	double dE = EE[0]-EE[3];
	double TT = DD[1]*v;
	// k1
	r11k1 = -2*TT*rho[1];
	r12k1 = rho12im*dE+TT*(2*rho[0]-1);
	i12k1 = -rho[1]*dE;
	// k2
	r11k2 = -2*TT*(rho[1]+r12k1*dT/2);
	r12k2 = (rho12im+i12k1*dT/2)*dE+TT*(2*(rho[0]+r11k1*dT/2)-1);
	i12k2 = -(rho[1]+r12k1*dT/2)*dE;
	// k3
	r11k3 = -2*TT*(rho[1]+r12k2*dT/2);
	r12k3 = (rho12im+i12k2*dT/2)*dE+TT*(2*(rho[0]+r11k2*dT/2)-1);
	i12k3 = -(rho[1]+r12k2*dT/2)*dE;
	// k4
	r11k4 = -2*TT*(rho[1]+r12k3*dT);
	r12k4 = (rho12im+i12k3*dT)*dE+TT*(2*(rho[0]+r11k3*dT)-1);
	i12k4 = -(rho[1]+r12k3*dT)*dE;
	//final
	rho[0] = rho[0] + (r11k1+2*r11k2+2*r11k3+r11k4)*dT/6;
	rho[1] = rho[1] + (r12k1+2*r12k2+2*r12k3+r12k4)*dT/6;
	rho12im = rho12im + (i12k1+2*i12k2+2*i12k3+i12k4)*dT/6;
	rho[2] = rho[1];
	rho[3] = 1-rho[0];

	rho_dot[0] = -2*TT*rho[1];
	rho_dot[3] = -rho_dot[0];
}

#endif
