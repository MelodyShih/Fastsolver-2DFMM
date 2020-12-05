#include <iostream>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include "boxtree.hpp"

#define rand01() ((double)rand()/RAND_MAX)

using namespace std;

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
	int p = (int) round(log(tol)/log(pow(2,0.5)/(4-pow(2,0.5))));
	return p;
}



int main(int argc, char *argv[])
{
	if (argc<2) {
		fprintf(stderr,
			"Usage: driver N tol\n"
			"Arguments:\n"
			"	N: number of charges\n"
		    "   tol: tolerance\n"
			"   totallevel: L\n"
			// "   p: p\n"
		);
		return 1;
	}
	int N;
	sscanf(argv[1],"%d",&N);
	double w, tol;
	sscanf(argv[2],"%lf",&w); tol=(double)w;
	int totallevel;
	sscanf(argv[3],"%d",&totallevel);
	// int p=4;
	// sscanf(argv[4],"%d",&p);
	int p = chooseP(tol);

	complex<double>* x = (complex<double>*) malloc(N * sizeof(complex<double>));
	double* q = (double*) malloc(N * sizeof(double));
	complex<double>* utrue  =(complex<double>*)
		                                   malloc(N * sizeof(complex<double>));
	complex<double>* uapprox=(complex<double>*)
		                                   malloc(N * sizeof(complex<double>));

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
#if 1
	// make N data
	// complex<double> array : x, y, q
	for(int i=0; i<N; i++){
		x[i] = complex<double>(rand01(), rand01());
		q[i] = 10*rand01();
		// q[i] = 1.0;
	}
	// exact(N, x, q, utrue);

	int rootlevel=0;
	int rootnum=0;

	Box* rootbox = new Box(rootlevel,0,0,rootnum,p);
	rootbox->buildtree(totallevel);
	cout<<"buildneighborinteractionlist"<<endl;
	rootbox->treetraverse(1);//buildneighborinteractionlist
	cout<<"buildTofo"<<endl;
	rootbox->treetraverse(4);//buildTofo
	cout<<"buildTifi"<<endl;
	rootbox->treetraverse(5);//buildTifi
	cout<<"buildTifo"<<endl;
	rootbox->treetraverse(6);//buildTifo
	cout<<"assign charges to box"<<endl;
	rootbox->assignchargestobox(totallevel, N, x);
	cout<<"upwardpass"<<endl;
	rootbox->upwardpass(totallevel, x, q);
	cout<<"downwardpass"<<endl;
	rootbox->downwardpass();
	//cout<<"printinteractionlist"<<endl;
	//rootbox->treetraverse(3);//buildneighborinteractionlist
	cout<<"build actual potential"<<endl;
	rootbox->buildactualpotential(totallevel, x, q, uapprox);

	// cout<<"exact potential:"<<endl;
	// printresult(N, utrue);
	// cout<<endl;
	cout<<"approximation:"<<endl;
	printresult(N, uapprox);

// // use the whole list
// 	double norm_true = norm_inf_vec(N, utrue);
// 	double norm_diff = norm_inf_diff(N, utrue, uapprox);
// 	// double norm_true = norm_two_vec(N, utrue);
// 	// double norm_diff = norm_two_diff(N, utrue, uapprox);

	// use one xi
	int idx = rand() % N;
	cout<<"idx: "<<idx<<endl;
	double norm_true = abs(exact_one(N, x, q, idx));
	double norm_diff = abs(exact_one(N, x, q, idx)-uapprox[idx].real());

	double rel_err = norm_diff/norm_true;

	cout<<"norm_true: "<<norm_true<<", norm_diff: "<<norm_diff<<endl;

	cout<<"p: "<<p<<endl;
	cout<<"relative error: "<<rel_err<<", tol="<<tol<<endl;
#endif
	return 0;
}
