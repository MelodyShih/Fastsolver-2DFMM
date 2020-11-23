#include <iostream> 
#include <cmath> 
using namespace std; 
class Box;

class Node
{
public:
    Box    *data=NULL;
    Node   *next=NULL;
    double *Tifo_dat=NULL;	
};

class Box
{
public:
	int i;// c_tau has coordinate i*h_level, h_level= 1/(pow(2,level))
	int j;
	int num;
	int level;

	double* q;
	double* x;
	double* y;

	double* qhat;
	double* uhat;

	double* Tofo_mat; // on parent, e.g. Tofo_mat is the map from
	                  // self (botleft) to its parent
	double* Tifi_mat; // on leaf, i.e. Tifi_mat is the map from 
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

	void computeoutcomingexp();
	void computeincomingexp();

	void buildTofo();
	void buildTifi();
	void buildTifo();
};
