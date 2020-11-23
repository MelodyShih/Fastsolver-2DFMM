#include <iostream>
#include "boxtree.hpp"

using namespace std;

/* Given a reference (pointer to pointer)
to the head of a list and an int,
inserts a new node  on the fr ont of the list. */
void push(Node** head_ref, Box* new_data)
{
    Node* new_node = new Node();
    new_node->data = new_data;
    new_node->next = (*head_ref);
    (*head_ref) = new_node;
}

void Box::buildtree(int level)
{
    if(level==0)
        return;
	int i = this->i;
	int j = this->j;
    this->topleft  = new Box(this->level+1, 2*i+1, 2*j, 2);
    this->topleft->parent = this;

    this->topright = new Box(this->level+1, 2*i+1, 2*j+1, 3);
    this->topright->parent = this;

    this->botleft  = new Box(this->level+1, 2*i, 2*j, 0);
    this->botleft->parent = this;

    this->botright = new Box(this->level+1, 2*i, 2*j+1, 1);
    this->botright->parent = this;

    this->topright->nextsibling = this->botleft;
    this->topleft->nextsibling  = this->topright;
    this->botright->nextsibling = this->topleft;
    this->botleft->nextsibling  = this->botright;

    this->topleft->buildtree(level-1);
    this->topright->buildtree(level-1);
    this->botleft->buildtree(level-1);
    this->botright->buildtree(level-1);
};

void Box::buildneighborinteractionlist()
{
	this->neighbor    = new Node();
	this->interaction = new Node();
	Box*  boxnow;
	Box*  boxparentnow;
	Box*  boxtoadd;
	boxnow = this;

	// Case 1: siblings share same parent (all must be in neighbor)
	for(int i=0; i<3; i++){
		if(boxnow->nextsibling == NULL){
			break;
		}
		push(&this->neighbor, boxnow->nextsibling);
		boxnow = boxnow->nextsibling; 
	}

	// Case 2: boxes from parent siblings (all must be in either neighbor or
	//         interaction)
	if(this->parent != NULL){
		boxnow = this->parent;
		for(int i=0; i<3; i++){
			if(boxnow->nextsibling == NULL){
				break;
			}

			//parent sibling's topleft
			boxtoadd = boxnow->nextsibling->topleft;
			if(abs(this->i-boxtoadd->i) <= 1 && abs(this->j-boxtoadd->j) <= 1){
				push(&this->neighbor, boxtoadd);
			}else{
				push(&this->interaction, boxtoadd);
			}

			//topright
			boxtoadd = boxnow->nextsibling->topright;
			if(abs(this->i-boxtoadd->i) <= 1 && abs(this->j-boxtoadd->j) <= 1){
				push(&this->neighbor, boxtoadd);
			}else{
				push(&this->interaction, boxtoadd);
			}

			//botleft
			boxtoadd = boxnow->nextsibling->botleft;
			if(abs(this->i-boxtoadd->i) <= 1 && abs(this->j-boxtoadd->j) <= 1){
				push(&this->neighbor, boxtoadd);
			}else{
				push(&this->interaction, boxtoadd);
			}

			//botright
			boxtoadd = boxnow->nextsibling->botright;
			if(abs(this->i-boxtoadd->i) <= 1 && abs(this->j-boxtoadd->j) <= 1){
				push(&this->neighbor, boxtoadd);
			}else{
				push(&this->interaction, boxtoadd);
			}

			boxnow = boxnow->nextsibling; 
		}

		// Case 3: rest of them
		if(this->parent->parent != NULL){
			boxparentnow = this->parent->parent;
			for(int i=0; i<3; i++){
				// if at root, break
				if(boxparentnow->nextsibling == NULL){
					break;
				}
				boxnow = boxparentnow->nextsibling->botleft; // go down one level
				for(int j=0; j<4; j++){
					// check that if boxnow's sibling touch this box's parent
					assert(boxnow->level == this->parent->level);
				    if(abs(this->parent->i - boxnow->i) <= 1 &&
				       abs(this->parent->j - boxnow->j) <= 1){
			            //parent sibling's topleft
			            boxtoadd = boxnow->topleft;
			            if(abs(this->i-boxtoadd->i) <= 1 && 
						   abs(this->j-boxtoadd->j) <= 1){
			            	push(&this->neighbor, boxtoadd);
			            }else{
			            	push(&this->interaction, boxtoadd);
			            }

			            //topright
			            boxtoadd = boxnow->topright;
			            if(abs(this->i-boxtoadd->i) <= 1 && 
						   abs(this->j-boxtoadd->j) <= 1){
			            	push(&this->neighbor, boxtoadd);
			            }else{
			            	push(&this->interaction, boxtoadd);
			            }

			            //botleft
			            boxtoadd = boxnow->botleft;
			            if(abs(this->i-boxtoadd->i) <= 1 && 
						   abs(this->j-boxtoadd->j) <= 1){
			            	push(&this->neighbor, boxtoadd);
			            }else{
			            	push(&this->interaction, boxtoadd);
			            }

			            //botright
			            boxtoadd = boxnow->botright;
			            if(abs(this->i-boxtoadd->i) <= 1 && 
						   abs(this->j-boxtoadd->j) <= 1){
			            	push(&this->neighbor, boxtoadd);
			            }else{
			            	push(&this->interaction, boxtoadd);
			            }
				    }
					boxnow = boxnow->nextsibling;
				}
				boxparentnow = boxparentnow->nextsibling; 
			}
		}
	}
}

void Box::printneighborlist()
{
	Node* now = this->neighbor;
	cout<<"neighbor of box ("<<this->i<<","<<this->j<<")"
		<<" on level "<<this->level
		<<endl;
	int n=0;
	while(now != NULL){
		if(now->data == NULL)
			break;
		for(int l=0; l<this->level; l++)
			cout<<"  ";
		cout<<"("<<now->data->i<<","<<now->data->j<<")"<<endl;
		now = now->next;
		n++;
	}
	cout<<"=>interaction of box ("<<this->i<<","<<this->j<<")"
		<<" on level "<<this->level
		<<", #interaction boxes "<<n
		<<endl;
}

void Box::printinteractionlist()
{
	Node* now = this->interaction;
	//if(this->level==3 && this->i == 2 && this->j == 3){
	int n=0;
	while(now != NULL){
		if(now->data == NULL)
			break;
		for(int l=0; l<this->level; l++)
			cout<<"  ";
		cout<<"("<<now->data->i<<","<<now->data->j<<")"<<endl;
		now = now->next;
		n++;
	}
	cout<<"=>interaction of box ("<<this->i<<","<<this->j<<")"
		<<" on level "<<this->level
		<<", #interaction boxes "<<n
		<<endl;
	//}
}

void performeaction(int action, Box* box)
{
	switch (action)
	{
    	case 1:
        	box->buildneighborinteractionlist();
			break;
    	case 2:
        	box->printneighborlist();
        	break;
		case 3:
        	box->printinteractionlist();
        	break;
    	default: // code to be executed if n doesn't match any cases
			cout<<"box ("<<box->i<<","<<box->j<<")"
				<<" on level "<<box->level
				<<endl;
	}
}

void Box::treetraverse(int action)
{
	if(this->botleft == NULL){ 
		performeaction(action, this);
		return;
	}
	this->botleft->treetraverse(action);
	this->botright->treetraverse(action);
	this->topleft->treetraverse(action);
	this->topright->treetraverse(action);
	performeaction(action, this);
}
