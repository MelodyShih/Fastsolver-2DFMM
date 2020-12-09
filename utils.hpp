#include <iostream>
#include <cmath>
#include <complex>

using namespace std;

void exact(int N, complex<double>* x, double* q, complex<double>* utrue);
double exact_one(int N, complex<double>* x, double* q, int i);
double norm_inf_vec(int N, complex<double>* utrue);
double norm_two_vec(int N, complex<double>* utrue);
double norm_inf_diff(int N, complex<double>* utrue, complex<double>* uapprox);
double norm_two_diff(int N, complex<double>* utrue, complex<double>* uapprox);
int chooseP(double tol);
complex<double> G(complex<double> x, complex<double> y);
void printresult(int N, complex<double>* u, complex<double>* uapprox);
