#include <iostream>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include "boxtree.hpp"

#define rand01() ((double)rand()/RAND_MAX)

using namespace std;

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

	// make N data
	// complex<double> array : x, y, q
	for(int i=0; i<N; i++){
		x[i] = complex<double>(rand01(), rand01());	
		q[i] = 10*rand01();
		cout << x[i] << endl;
	}

	int totallevel=3; //L, 0,1,2, ..., L-1, determine by number of charges
	int rootlevel=0;
	int rootnum=0;
	Box* rootbox = new Box(rootlevel,0,0,rootnum,p);
	rootbox->buildtree(totallevel); 
	rootbox->assignchargestobox(totallevel, N, x);

	//rootbox->treetraverse(1);//build neighborinterationlist
	//rootbox->treetraverse(2);//print neighborlist
	//rootbox->treetraverse(3);//print interactionlist
	//rootbox->treetraverse(4);//build Tofo

	//cout<<"Downward pass:"<<endl;
	//rootbox->downwardpass(0);
	//cout<<endl;
	//cout<<"Upward pass:"<<endl;
	//rootbox->upwardpass(0);

	return 0;
}
