#include <iostream> 
#include <cmath> 
using namespace std; 
class Box;

class Node
{
	public:
		Box    *data=NULL;
		Node   *next=NULL;
		double *Tifo_dat=NULL;// pxp matrix
};

class Box
{
	public:
		int i;// xcoo
		      // c_tau has coordinate (i+0.5)*h_level, h_level= 1/(pow(2,level))
		int j;// ycoo

		double cx;
		double cy;

		int num;
		int level;

		int* idxcharge; //q[idxcharge] mpoints
		//double* q;
		//double* x;
		//double* y;

		double* qhat; // px1 vector
		double* uhat; // px1 vector

		double* Tofo_mat; // pxp matrix, on leaf, e.g. Tofo_mat is the map from
		// self (botleft) to its parent
		double* Tifi_mat; // pxp matrix, on leaf, i.e. Tifi_mat is the map from 
		                  // parent to self	

		Box* parent;
		Box* topleft;
		Box* topright;
		Box* botleft;
		Box* botright;

		Box*  nextsibling;
		Node* neighbor;
		Node* interaction;

		Box(int level, int i, int j, int num)
		{
			this->i = i;
			this->j = j;
			this->num = num;
			this->level = level;
			this->parent   = NULL;
			this->topleft  = NULL; 
			this->topright = NULL; 
			this->botleft  = NULL; 
			this->botright = NULL; 
			this->nextsibling = NULL;
			//cout<<"(level, num)=("<<level<<","<<i<<","<<j<<")"<<endl;
		};

		void buildtree(int numlevel);
		void buildneighborinteractionlist();
		void printneighborlist();
		void printinteractionlist();
		void treetraverse(int action);
		void downwardpass(int action);
		void upwardpass(int action);
		void buildactualpotential();

		void computeoutcomingexp();//TODO: rename
		void computeincomingexp();

		void buildTofo();
		void buildTifi();
		void buildTifo();
};
