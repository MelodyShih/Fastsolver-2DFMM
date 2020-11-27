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

int main(int argc, char *argv[])
{
	if (argc<2) {
		fprintf(stderr,
			"Usage: driver N tol\n"
			"Arguments:\n"
			"	N: number of charges\n"
		    "   tol: tolerance\n");
		return 1;
	}
	int N;
	sscanf(argv[1],"%d",&N);
	double w, tol;
	sscanf(argv[6],"%lf",&w); tol=(double)w;
	int p = 4;

	complex<double>* x = (complex<double>*) malloc(N * sizeof(complex<double>)); 
	double* q = (double*) malloc(N * sizeof(double)); 
	complex<double>* utrue  =(complex<double>*) 
		                                   malloc(N * sizeof(complex<double>));
	complex<double>* uapprox=(complex<double>*) 
		                                   malloc(N * sizeof(complex<double>));
	
	// make N data
	// complex<double> array : x, y, q
	for(int i=0; i<N; i++){
		x[i] = complex<double>(rand01(), rand01());	
		q[i] = 10*rand01();
	}
	exact(N, x, q, utrue);

	int totallevel=3; //L, 0,1,2, ..., L-1, determine by number of charges
	int rootlevel=0;
	int rootnum=0;

	Box* rootbox = new Box(rootlevel,0,0,rootnum,p);
	rootbox->buildtree(totallevel); 
	rootbox->treetraverse(1);//buildneighborinteractionlist
	rootbox->treetraverse(4);//buildTofo
	rootbox->treetraverse(5);//buildTifi
	rootbox->treetraverse(6);//buildTifo
	rootbox->assignchargestobox(totallevel, N, x);
	rootbox->upwardpass(totallevel, x, q);
	rootbox->downwardpass();
	rootbox->buildactualpotential(totallevel, x, q, uapprox);

	cout<<"exact potential:"<<endl;
	printresult(N, utrue);
	cout<<endl;
	cout<<"approximation:"<<endl;
	printresult(N, uapprox);
	return 0;
}
