#include "mat2.h"
#include "pot.h"
#include "utilities.h"
#include <iostream>

using namespace std;

int main()
{
	mat2 (*Hp)(double)=&H1;
	mat2 vec,val;
	for (double x=-10; x<=10; x+=0.1)
	{
		Hp(x).eig(vec,val);
		cout<<x<<'\t'<<val[0]<<'\t'<<val[3]<<endl;
	}
	return 0;
}
	
