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
		    "   tol: tolerance\n"
			"   totallevel: L\n");
		return 1;
	}
	int N;
	sscanf(argv[1],"%d",&N);
	double w, tol;
	sscanf(argv[2],"%lf",&w); tol=(double)w;
	int p = 4;
	int totallevel;
	sscanf(argv[3],"%d",&totallevel);

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
	cout<<"build actual potential"<<endl;
	rootbox->buildactualpotential(totallevel, x, q, uapprox);

	cout<<"exact potential:"<<endl;
	printresult(N, utrue);
	cout<<endl;
	cout<<"approximation:"<<endl;
	printresult(N, uapprox);
	return 0;
}
