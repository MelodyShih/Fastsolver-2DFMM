#include <iostream> 
#include <cmath> 
using namespace std; 
class Box;

class Node
{
public:
    Box  *data=NULL;
    Node *next=NULL;
};

class Box
{
public:
	int i;
	int j;
	int num;
	int level;
	double* q;

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

	void buildtree(int level);
	void buildneighborinteractionlist();
	void printneighborlist();
	void printinteractionlist();
	void treetraverse(int action);
};
