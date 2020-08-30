#include "mat2.h"
#include "pot.h"
#include "utilities.h"
#include <iostream>
#include <cmath>
#include <time.h>
#include <mpi.h>

using namespace std;

int main()
{
	MPI_Init(NULL,NULL);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	//============
	srand(time(0)+rank*10);
	mat2 (*Hp)(double) = &H2;
	double m=2000.0;
	double x0 = -10;
	double dTe = 1/0.05/80;
	int num_iter = 600;
	int num_k = 60;
	double E_min = exp(-4), E_max = exp(1);
	mat2 rho0;

	rho0[0]=1;rho0[1]=0;rho0[2]=0;rho0[3]=0;

	for (int tk=0; tk<num_k; tk++)
	{
		double kk=sqrt(2*m*E_min*pow(E_max/E_min,tk/(double)(num_k-1)));
		double v0 = kk/m;
		double dTa = 1/v0/40;
		int atom_time_scale = (int)(dTa/dTe);
		if (atom_time_scale == 0)
		{
			cout<<"Error, atomic time scale too small"<<endl;
			exit(EXIT_FAILURE);
		}
		double trans1=0,trans2=0,reflect1=0;
		double btrans1=0,btrans2=0,breflect1=0;
		for (int t1=0; t1<num_iter; t1++)
		{
			double xx = x0 + (rand()/(double)RAND_MAX)*2-1;
			double vv = v0;
			int s=0;
			mat2 Dir,Ene,Vec;
			mat2 rho = rho0;
			mat2 rho_dot;
			double rho12im = 0;
			int count = 0;
			double ek_res;
			while(fabs(xx) < fabs(x0)+1.1)
			{
				// do atomic step
				if (count%atom_time_scale == 0)
				{
					x_v_rk4(dTa,xx,vv,m,s,Hp);
					Hp(xx).eig(Vec,Ene);
					Dir = DD(xx,Hp);
				}
				// electronic step
				rho_rk4(dTe,rho,rho12im,rho_dot,vv,Dir,Ene);
				ek_res = m*vv*vv/2+Ene[s]-Ene[3-s];
				if (ek_res > 0 && -rho_dot[s]*dTe/rho[s] > rand()/(double)RAND_MAX)
				{
					if (vv>=0)
						vv = sqrt(2*ek_res/m);
					else
						vv = -sqrt(2*ek_res/m);
					s = 3-s;
				}
				count++;
			}
			if(xx < 0 && s == 0)
				reflect1++;
			else if(s == 0)
				trans1++;
			else
				trans2++;
		}
		MPI_Allreduce(&trans1,&btrans1,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		MPI_Allreduce(&trans2,&btrans2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		MPI_Allreduce(&reflect1,&breflect1,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		if (rank == 0)
			cout<<log(kk*kk/2/m)<<'\t'<<btrans1/num_iter/size<<'\t'<<btrans2/num_iter/size<<'\t'<<breflect1/num_iter/size<<endl;
	}

	MPI_Finalize();
	return 0;
}
	
