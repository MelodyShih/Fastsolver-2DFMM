#include <iostream>
#include <stdlib.h>     /* srand, rand */
#include <chrono>
#include "utils.hpp"

#define rand01() ((double)rand()/RAND_MAX)

using namespace std;

complex<double> G(complex<double> x, complex<double> y)
{
	return log(x-y);
}

void printresult(int N, complex<double>* u, complex<double>* uapprox){
	for(int i=0; i<N; i++){
		cout<<u[i]<<" "<<uapprox[i]<<endl;
	}
}

void exact(int N, complex<double>* x, double* q, complex<double>* utrue)
{
	for(int i=0; i<N; i++){
		utrue[i] = complex<double>(0,0);
		for(int j=0; j<N; j++){
			if(i == j) continue;
			utrue[i] = utrue[i] + G(x[i], x[j])*q[j];
		}
	}
}

double exact_one(int N, complex<double>* x, double* q, int i)
{// exact real part of certain utrue[i]
		complex<double> utruei = complex<double>(0,0);
		for(int j=0; j<N; j++){
			if(i == j) continue;
			utruei = utruei + G(x[i], x[j])*q[j];
		}
		// cout<<"true u[idx]:"<<utruei<<endl;
		double utrueireal = utruei.real();
	return utrueireal;
}

double norm_inf_vec(int N, complex<double>* utrue)
{ // compute the L_inf norm
	double norm_inf = 0;
	for(int i=0; i<N; i++){
			if(abs(utrue[i].real()) > norm_inf){
				norm_inf = abs(utrue[i].real());
			}
	}
	return norm_inf;
}

double norm_two_vec(int N, complex<double>* utrue)
{ // compute the L2 norm
	double norm2 = 0;
	for(int i=0; i<N; i++){
		norm2 += pow(utrue[i].real(), 2);
	}
	return pow(norm2, 0.5);
}

double norm_inf_diff(int N, complex<double>* utrue, complex<double>* uapprox)
{ // compute the L_inf norm of difference between two vectors
	double norm_inf = 0;
	for(int i=0; i<N; i++){
			if(abs(utrue[i].real()-uapprox[i].real()) > norm_inf){
				norm_inf = abs(utrue[i].real()-uapprox[i].real());
			}
	}
	return norm_inf;
}

double norm_two_diff(int N, complex<double>* utrue, complex<double>* uapprox)
{ // compute the L2 norm of difference between two complex vectors
	double norm2 = 0;
	for(int i=0; i<N; i++){
		norm2 += pow(utrue[i].real()-uapprox[i].real(), 2);
	}
	return pow(norm2,0.5);
}

int chooseP(double tol)
{ // compute length P given tolerance
	return (int) round(log(tol)/log(pow(2,0.5)/(4-pow(2,0.5))));
}

#if 0
	// testing expansion of log
	complex<double> xx(1.0, 0.0);
	complex<double> yy(0.0, 0.2);
	complex<double> cc(0.1, 0.1);
	complex<double> approx(0.0, 0.0);
	cout<<"exact:" << G(xx,yy)<<endl;
	approx = G(xx,cc);
	for(int j=1; j<p; j++){
		approx += ( -1./ (double) j)*pow(yy-cc, j)*pow(xx-cc,-j);
	}
	cout<<"approx:"<<approx<<endl;
#endif
