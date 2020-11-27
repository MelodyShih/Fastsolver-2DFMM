#include <iostream>
#include <cmath>
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

void Box::buildtree(int numlevel)
{
	if(numlevel==1)
		return;
	int i = this->i;
	int j = this->j;
	int level = this->level;
	int p = this->p;
	double h = 1/pow(2, level);

	this-> c = complex<double>((i+0.5)*h, (j+0.5)*h);

	this->topleft  = new Box(this->level+1, 2*i+1, 2*j, 2, p);
	this->topleft->parent = this;

	this->topright = new Box(this->level+1, 2*i+1, 2*j+1, 3, p);
	this->topright->parent = this;

	this->botleft  = new Box(this->level+1, 2*i, 2*j, 0, p);
	this->botleft->parent = this;

	this->botright = new Box(this->level+1, 2*i, 2*j+1, 1, p);
	this->botright->parent = this;

	this->topright->nextsibling = this->botleft;
	this->topleft->nextsibling  = this->topright;
	this->botright->nextsibling = this->topleft;
	this->botleft->nextsibling  = this->botright;

	this->topleft->buildtree(numlevel-1);
	this->topright->buildtree(numlevel-1);
	this->botleft->buildtree(numlevel-1);
	this->botright->buildtree(numlevel-1);
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

void Box::assignchargestobox(int totallevel, int N, complex<double>* x)
{
	int  numleafbox = (int) pow(4,totallevel-1);
    int* idxcharge[numleafbox];
    int* numchargeperleafbox = (int*) malloc(numleafbox*sizeof(int));
    int* idxleafbox = (int*) malloc(numleafbox*sizeof(N));
    double boxsize = 1/pow(2,totallevel-1);
	cout <<"boxsize "<<boxsize<<endl;

    for(int b=0; b<numleafbox; b++){
        numchargeperleafbox[b] = 0.0;
    }

    for(int i=0; i<N; i++){
        double cx = x[i].real();
        double cy = x[i].imag();
        int ix = floor(cx/boxsize);
        int iy = floor(cy/boxsize);
        int idxbox = iy*pow(2,totallevel-1) + ix;

        idxleafbox[i] = numchargeperleafbox[idxbox];
        numchargeperleafbox[idxbox] += 1;
    }

    for(int b=0; b<numleafbox; b++){
        idxcharge[b] = (int*) malloc((numchargeperleafbox[b]+1)*sizeof(int));
        idxcharge[b][numchargeperleafbox[b]] = -1;
    }

    for(int i=0; i<N; i++){
        double cx = x[i].real();
        double cy = x[i].imag();
        int ix = floor(cx/boxsize);
        int iy = floor(cy/boxsize);
        int idxbox = iy*pow(2,totallevel-1) + ix;
        idxcharge[idxbox][idxleafbox[i]] = i;
    }
#if 0
    for(int b=0; b<numleafbox; b++){
		for(int i=0; i<numchargeperleafbox[b]+1; i++){
			cout << idxcharge[b][i]<<" ";
		}
		cout<<endl;
    }
#endif
	this->assignidxtoleaf(totallevel-1, idxcharge);
}

void performaction(int action, Box* box)
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
		case 4:
			box->buildTofo();	
			break;
		case 5:
			box->buildTifi();	
			break;
		case 6:
			box->buildTifo();	
			break;
		default: // code to be executed if n doesn't match any cases
			for(int l=0; l<box->level; l++)
				cout<<"  ";
			cout<<"box ("<<box->i<<","<<box->j<<")"<<endl;
	}
}

void Box::treetraverse(int action)
{
	if(this->botleft == NULL){ 
		performaction(action, this);
		return;
	}
	this->botleft->treetraverse(action);
	this->botright->treetraverse(action);
	this->topleft->treetraverse(action);
	this->topright->treetraverse(action);
	performaction(action, this);
}

void Box::assignidxtoleaf(int level, int** idxchargearray)
{
	if(this->level == level){
		int idxbox = this->i + this->j*(int)pow(2,level);
		this->idxcharge = idxchargearray[idxbox];
		if(idxchargearray[idxbox][0] == -1) 
			return;
		i = 0;
		while(idxchargearray[idxbox][i] != -1){
			cout << this->idxcharge[i] << " ";
			i++;
		}
		cout << endl;
		return;
	}
	this->botleft ->assignidxtoleaf(level, idxchargearray);
	this->botright->assignidxtoleaf(level, idxchargearray);
	this->topleft ->assignidxtoleaf(level, idxchargearray);
	this->topright->assignidxtoleaf(level, idxchargearray);
}

void Box::downwardpass(int action)
{
	if(this->botleft == NULL){ 
		return;
	}
	performaction(action, this); //TODO this->computeincomingexp()
	this->botleft ->downwardpass(action);
	this->botright->downwardpass(action);
	this->topleft ->downwardpass(action);
	this->topright->downwardpass(action);
}

void Box::upwardpass(int action)
{
	if(this->botleft == NULL){ 
		performaction(action, this); // TODO  this->computeoutgoingexp();
		return;
	}
	this->botleft ->upwardpass(action);
	this->botright->upwardpass(action);
	this->topleft ->upwardpass(action);
	this->topright->upwardpass(action);
	performaction(action, this);
}

void Box::computeoutgoingexp()
{
	// Case 1: Box is on the leaf of the tree (this->level == L-1)
	// Apply outgoing-from-sources map T_tau^{ofs}, see (7.2)


	// Case 2: Box is a parent of 4 nodes
	// Apply outgoing from outgoing map T_tau^{ofo}, see (7.3)
	// (get from looping through its child)
	// this->topleft->T_tau^{ofo} x this->topleft->qhat
}

void Box::computeincomingexp()
{
	// Case 1: level 0, 1
	// uhat = 0


	// Case 2: all other levels (7.4) (7.5)
	// Apply T_sigma,parent^{ifi} to uhat_parent


	// Loop through all the boxes in the interaction list sigma
	// Apply T_sigma,interaction^{ifo} to qhat_interaction
}

void Box::buildTifi()
{
	//Theorem 7.2
	//c_sigma (this): (this->cx, this->cy)
	//c_tau (this->parent): (this->parent->cx, this->parent->cy)
}

int fact(int n)
{
	// Returns factorial of n
    int res = 1;
    for (int i = 2; i <= n; i++)
        res = res * i;
    return res;
}
int nCr(int n, int r)
{
    return fact(n) / (fact(r) * fact(n - r));
}

void Box::buildTofo()
{
	if(this->level == 0) return; // do nothing for root box
	//Theorem 7.1
	complex<double> cparent = this->parent->c;
	complex<double> c = this->c;

	int p = this->p;
	this->Tofo_mat = (complex<double>*) malloc(p*p * sizeof(complex<double>));
	cout<<"Tofo for box " <<this->i<<","<<this->j<<" on level "<<this->level<<endl;
	for(int j=0; j<p; j++){
		for(int i=0; i<p; i++){
			int idx = j*p+i;
			if (i<=j){
				this->Tofo_mat[idx] = (double) nCr(j,i)*
					                        pow(cparent-c, j-i); 
			}else{
				this->Tofo_mat[idx] = complex<double>(0,0); 
			}
			cout<<this->Tofo_mat[idx]<<" ";
		}
		cout<<endl;
	}
}

void Box::buildTifo()
{
	//Theorem 6.26 
	//Loop through all boxes in this->interaction
	//now = this->interaction
	//now->data->cx, now->data->cy
	//now = now->next
	//until now == NULL
}
