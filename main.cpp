/*
0.0001
k =1
n=1
omega jako pierwsza rezonansowa
energia powinna wyniosić 4,93
N = 100
w = 2 pi^2 / ...
*/

#include <cstring>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <fstream>
#include <cstdlib>
#define _USE_MATH_DEFINES
using namespace std;
const int N = 100;
const double t=0.0001;
const double w = 3.0/2.0*M_PI*M_PI*0;
const int k =1;

void H_t(double T,double Y[],double H[],double x[])
{
	H[0]=0;
	H[101]=0;
	for(int ii=1;ii<N;ii++)
	{
		H[ii]=(double)-N*N/2*(Y[ii-1]+Y[ii+1]-2*Y[ii])+k*(x[ii]-1.0/2)*Y[ii]*sin(w*T);
	}
}

void step(double Yr[],double Yi[],double Hr[],double Hi[],double x[], int num)
{
	H_t(num*t,Yi,Hi,x);
	for(int ii=0;ii<N+1;ii++) Yr[ii]=Yr[ii]+Hi[ii]*t/2;
	H_t(num*t+t/2,Yr,Hr,x);
	for(int ii=0;ii<N+1;ii++) Yi[ii]=Yi[ii]-Hr[ii]*t;
	H_t((num+1)*t,Yi,Hi,x);
	for(int ii=0;ii<N+1;ii++) Yr[ii]=Yr[ii]+Hi[ii]*t/2;
}

double E(double Yr[],double Yi[],double Hr[],double Hi[])
{
	double tmp;
	double sum=0;
	for(int ii=0;ii<N+1;ii++)
	{
		tmp=Hr[ii]*Yr[ii]+Hi[ii]*Yi[ii];
		sum+=tmp/N;
	}
	return sum;
}

double pozX(double Yr[],double Yi[],double Hr[],double Hi[],double x[])
{
	double tmp;
	double sum=0;
	for(int ii=0;ii<N+1;ii++)
	{
		tmp=x[ii]*(Yr[ii]*Yr[ii]+Yi[ii]*Yi[ii]);
		sum+=tmp/N;
	}
	return sum;
}

double Norm(double Yr[],double Yi[],double Hr[],double Hi[])
{
	double tmp;
	double sum=0;
	for(int ii=0;ii<N+1;ii++)
	{
		tmp=(Yr[ii]*Yr[ii]+Yi[ii]*Yi[ii]);
		sum+=tmp/N;
	}
	return sum;
}

int main()
{

	double x[N+1];
	double Yr[N+1];
	double Yi[N+1];
	int n=1;
	double Hr[N+1];
	double Hi[N+1];


	for(int ii=0;ii<101;ii++)
	{
		x[ii]=(double)ii/N;
		Yr[ii]=sqrt(2)*sin(n*M_PI*x[ii]);
		Yi[ii]=0;
		Hi[ii]=0;
	}

	H_t(0,Yr,Hr,x);

	double En = E(Yr,Yi,Hr,Hi);
	cout<<En<<endl;

	fstream dataOut;
	dataOut.open("data.txt", ios::out | ios::trunc);
	//dataOut<<"E\tx\tN"<<endl;

	int steps = 500000;
	for (int ii=0;ii<steps;ii++)
		{
			step(Yr,Yi,Hr,Hi,x,ii);
			if(ii%2==0)
			{
				dataOut<<ii<<"\t"<<E(Yr,Yi,Hr,Hi)<<"\t"<<pozX(Yr,Yi,Hr,Hi,x)<<"\t"<<Norm(Yr,Yi,Hr,Hi)<<endl;
			}
			if(ii%(steps/100)==0) cout<<ii*100/steps<<"%"<<endl;
		}


	cout<<"E= "<<E(Yr,Yi,Hr,Hi);
	cout<<"  x= "<<pozX(Yr,Yi,Hr,Hi,x);
	cout<<"  N= "<<Norm(Yr,Yi,Hr,Hi)<<endl;

	dataOut.close();


/*
* TODO zależność energii od czasu i krzywa rezonansowa z fitem
*/

	return 0;
}
