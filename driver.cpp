#include <iostream>
#include "boxtree.hpp"
using namespace std;

int main()
{
	int level=3;
	Box* rootbox = new Box(0,0,0,0);
    rootbox->buildtree(level); 

	rootbox->treetraverse(1);//build neighborinterationlist
	//rootbox->treetraverse(2);//print neighborlist
	rootbox->treetraverse(3);//print interactionlist

    return 0;
}