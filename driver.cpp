#include <iostream>
#include <stdlib.h>     /* srand, rand */
#include <chrono>
#include "boxtree.hpp"
#include "utils.hpp"

#define rand01() ((double)rand()/RAND_MAX)

using namespace std;

int main(int argc, char *argv[])
{
	if (argc<2) {
		fprintf(stderr,
			"Usage: driver N tol [p [totallevel]]\n"
			"Arguments:\n"
			"	N: number of charges\n"
		    "   tol: tolerance\n"
			"   totallevel: L\n"
			"   p: p\n"
		);
		return 1;
	}
	int N;
	sscanf(argv[1],"%d",&N);
	double w, tol;
	sscanf(argv[2],"%lf",&w); tol=(double)w;

	int p = chooseP(tol);
    if(argc>3){
	    sscanf(argv[3],"%d",&p);
    }

	int totallevel = ceil(log10((double)N/p)/log10(4.0))+1;
    //cout << pow(4,totallevel-1)*p <<endl;
    if (argc>4){
	    sscanf(argv[4],"%d",&totallevel);
    }
    cout << "[info] p = "<<p<<", totallevel = "<<totallevel<<endl;;

	double* q = (double*) malloc(N * sizeof(double));
	complex<double>* x = (complex<double>*) malloc(N * sizeof(complex<double>));
	complex<double>* utrue  =
                        (complex<double>*) malloc(N * sizeof(complex<double>));
	complex<double>* uapprox=
                        (complex<double>*) malloc(N * sizeof(complex<double>));

	// Create data
	for(int i=0; i<N; i++){
        uapprox[i] = complex<double>(0.0, 0.0);
		x[i] = complex<double>(rand01(), rand01());
		//q[i] = 0.01*rand01();
		q[i] = 1.0;
        //cout << x[i] <<endl;
	}
    q[N-1]=0.0;

	int rootlevel=0;
	int rootnum=0;

    auto start = chrono::steady_clock::now();
    double time;
    double precomptime=0.0;
	Box* rootbox = new Box(rootlevel,0,0,rootnum,p);
	rootbox->buildtree(totallevel);
    auto end = chrono::steady_clock::now();
    time = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    precomptime += time;
    cout <<"[time] Build tree            : "<<time<< " ms" << endl;
    start = chrono::steady_clock::now();
	rootbox->treetraverse(1);// buildneighborinteractionlist
    end = chrono::steady_clock::now();
    time = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    precomptime += time;
    cout <<"[time] Build nei, inter list : "<<time<< " ms" << endl;
    start = chrono::steady_clock::now();
	rootbox->treetraverse(4);// buildTofo, Tifi, Tifo
    end = chrono::steady_clock::now();
    time = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    precomptime += time;
    cout <<"[time] Build transfer op     : "<<time<< " ms" <<endl;
    cout <<"[time] Cost indep. of charges: "<<precomptime<< " ms" << endl;

    start = chrono::steady_clock::now();
	rootbox->assignchargestobox(totallevel, N, x);
    end = chrono::steady_clock::now();
    cout << "[time] Assign charges to box : "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;

    start = chrono::steady_clock::now();
	rootbox->upwardpass(totallevel, x, q);
	rootbox->downwardpass();
	rootbox->buildactualpotential(totallevel, x, q, uapprox);
    end = chrono::steady_clock::now();
    cout << "[time] Perform FMM           : "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;
    
    if(N<10000){
        start = chrono::steady_clock::now();
        exact(N, x, q, utrue);
        end = chrono::steady_clock::now();
        cout << "[time] Direct                : "
             << chrono::duration_cast<chrono::milliseconds>(end - start).count()
             << " ms" << endl;

        // use the whole list
        double norminf_true = norm_inf_vec(N, utrue);
        double norminf_diff = norm_inf_diff(N, utrue, uapprox);
        double norml2_true  = norm_two_vec(N, uapprox);
        double norml2_diff  = norm_two_diff(N, utrue, uapprox);
        cout<<"[acc ] relative l2   error   : "<<norml2_diff/norml2_true<<endl;
        cout<<"[acc ] relative linf error   : "<<norminf_diff/norminf_true<<endl;
    }

	// use one xi
	int idx = rand() % N;
	double norm_true = abs(exact_one(N, x, q, idx));
	double norm_diff = abs(exact_one(N, x, q, idx)-uapprox[idx].real());
	double rel_err = norm_diff/norm_true;

	cout<<"[acc ] norm_true: "<<norm_true<<", norm_diff: "<<norm_diff<<endl;
	cout<<"[acc ] relative error (1 idx): "<<rel_err<<" (idx="<<idx<<")"<<endl;
	return 0;
}
