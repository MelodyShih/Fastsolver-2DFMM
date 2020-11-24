#include <iostream>
#include "boxtree.hpp"
using namespace std;

int main()
{
	// make N data
	// double array : x, y, q
	int level=3; //L, 0,1,2, ..., L-1
	Box* rootbox = new Box(0,0,0,0);
	rootbox->buildtree(level); 

	rootbox->treetraverse(1);//build neighborinterationlist
	rootbox->treetraverse(2);//print neighborlist
	rootbox->treetraverse(3);//print interactionlist

	cout<<"Downward pass:"<<endl;
	rootbox->downwardpass(0);
	cout<<endl;
	cout<<"Upward pass:"<<endl;
	rootbox->upwardpass(0);

	return 0;
}
