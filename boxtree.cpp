#include <iostream>
#include <cmath>
#include <assert.h>
#include "boxtree.hpp"
#include "utils.hpp"

using namespace std;

double fact(int n)
{
    // Returns factorial of n
    double res = 1;
    for (int i = 2; i <= n; i++)
        res = res * i;
    return res;
}

double nCr(int n, int r)
{
    double result=1.0;
    for(int i=1; i<r; i++){
        result = result*(n-i+1)/(double) i;
    }
    //return result;
    return fact(n) / (fact(r) * fact(n - r));
}



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
            box->buildTifi();
            box->buildTifo();
            break;
        default: // code to be executed if n doesn't match any cases
            for(int l=0; l<box->level; l++)
                cout<<"  ";
            cout<<"box ("<<box->i<<","<<box->j<<")"<<endl;
    }
}

void Box::buildtree(int numlevel)
{
    if(numlevel==1)
        return;
    int i = this->i;
    int j = this->j;
    int level = this->level;
    int p = this->p;

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
    if(this->level == 0){
        return;
    }
    Box* boxnow = this;;
    for(int i=0; i<3; i++){
        push(&this->neighbor, boxnow->nextsibling);
        boxnow = boxnow->nextsibling;
    }
    if(this->level > 1){
        Node* now = this->parent->neighbor;
        Box* nowparent = now->data;
        Box* boxtoadd;
        while(nowparent != NULL){
            boxtoadd = nowparent->topleft;
            for(int i=0; i<4; i++){
                if(abs(this->i-boxtoadd->i) <= 1 && abs(this->j-boxtoadd->j) <= 1){
                    push(&this->neighbor, boxtoadd);
                }else{
                    push(&this->interaction, boxtoadd);
                }
                boxtoadd = boxtoadd->nextsibling;
            }
            now = now->next;
            nowparent = now->data;
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
    cout<<"=>neighbor of box ("<<this->i<<","<<this->j<<")"
        <<" on level "<<this->level
        <<", #neighbor boxes "<<n
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
    int* idxleafbox = (int*) malloc(N*sizeof(int));
    double boxsize = 1.0/pow(2,totallevel-1);
    //cout << boxsize <<endl;

    for(int b=0; b<numleafbox; b++){
        numchargeperleafbox[b] = 0.0;
    }

    for(int i=0; i<N; i++){
        double cx = x[i].real();
        double cy = x[i].imag();
        int ix = floor(cx/boxsize);
        int iy = floor(cy/boxsize);
        int idxbox = (int) (iy*pow(2,totallevel-1) + ix);

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
        int idxbox = (int) (iy*pow(2,totallevel-1) + ix);
        idxcharge[idxbox][idxleafbox[i]] = i;
    }
#if 0
    int n = 0;
    for(int b=0; b<numleafbox; b++){
        if( b % (int)pow(2,totallevel-1) == 0)
            cout<<endl;
        n += numchargeperleafbox[b];
        cout<<numchargeperleafbox[b]<<" ";
        //cout << "b="<<b<<endl;
        //for(int i=0; i<numchargeperleafbox[b]; i++){
        //	cout << idxcharge[b][i]<<" ";
        //}
        //cout<<endl;
    }
    cout << n <<endl;
#endif
    this->assignidxtoleaf(totallevel-1, idxcharge, numchargeperleafbox);
}

#if 0
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
#endif
void Box::treetraverse(int action)
{
    performaction(action, this);
    if(this->botleft == NULL){
        return;
    }
    this->botleft->treetraverse(action);
    this->botright->treetraverse(action);
    this->topleft->treetraverse(action);
    this->topright->treetraverse(action);
}

void Box::assignidxtoleaf(int level, int** idxchargearray, int* numchargeperleafbox)
{
    if(this->level == level){
        int idxbox = this->i + this->j*(int)pow(2,level);
        this->idxcharge = idxchargearray[idxbox];
        this->ncharge = numchargeperleafbox[idxbox];
#ifdef DEBUG
        cout<<this->i<<","<<this->j<<","<<this->ncharge<<endl;
        if(idxchargearray[idxbox][0] == -1)
            return;
        i = 0;
        while(idxchargearray[idxbox][i] != -1){
            cout << this->idxcharge[i] << " ";
            i++;
        }
        cout << endl;
#endif
        return;
    }
    this->botleft ->assignidxtoleaf(level, idxchargearray, numchargeperleafbox);
    this->botright->assignidxtoleaf(level, idxchargearray, numchargeperleafbox);
    this->topleft ->assignidxtoleaf(level, idxchargearray, numchargeperleafbox);
    this->topright->assignidxtoleaf(level, idxchargearray, numchargeperleafbox);
}

void Box::downwardpass()
{
    if(this->botleft == NULL){
        this->computeincomingexp();
        return;
    }
    this->computeincomingexp();
    this->botleft ->downwardpass();
    this->botright->downwardpass();
    this->topleft ->downwardpass();
    this->topright->downwardpass();
}

void Box::upwardpass(int totallevel, complex<double>* x, double* q)
{
    if(this->botleft == NULL){
        //performaction(0, this);
        this->computeoutgoingexp(totallevel, x, q);
        return;
    }
    this->botleft ->upwardpass(totallevel, x, q);
    this->botright->upwardpass(totallevel, x, q);
    this->topleft ->upwardpass(totallevel, x, q);
    this->topright->upwardpass(totallevel, x, q);
    //performaction(0, this);
    this->computeoutgoingexp(totallevel, x, q);
}

void Box::computeoutgoingexp(int totallevel, complex<double>* x, double* q)
{
    int p = this->p;
    this->qhat = (complex<double>*) malloc(p * sizeof(complex<double>));
    if(this->level == totallevel-1){
        // Case 1: Box is on the leaf of the tree (this->level == L-1)
        // Apply outgoing-from-sources map T_tau^{ofs}, see (7.2)
        for(int j = 0; j<p; j++){ // row
            this->qhat[j] = complex<double>(0,0);
            for(int i = 0; i < this->ncharge; i++){// for charges in this box
                if(j==0){// first element
                    this->qhat[j] += q[this->idxcharge[i]]; // how to pass q?
                }else{ // othere elements
                    this->qhat[j] += ( -1./ (double) j)
                        *pow(x[this->idxcharge[i]]-this->c, j)
                        *q[this->idxcharge[i]] ;
                }
            }
        }
    }else{
        // Case 2: Box is a parent of 4 nodes
        // Apply outgoing from outgoing map T_tau^{ofo}, see (7.3)
        // (get from looping through its child)
        // this->topleft->T_tau^{ofo} x this->topleft->qhat
        for(int j=0; j<p;j++){//row
            this->qhat[j] = complex<double>(0,0);
            for(int i=0; i<p; i++){//column
                int idx = j*p+i;
                double scale = (i==0? 1: -i);
                this->qhat[j] +=
                    this->topleft->Tofo_mat[idx]*this->topleft->qhat[i]*scale;
                this->qhat[j] +=
                    this->topright->Tofo_mat[idx]*this->topright->qhat[i]*scale;
                this->qhat[j] +=
                    this->botleft->Tofo_mat[idx]*this->botleft->qhat[i]*scale;
                this->qhat[j] +=
                    this->botright->Tofo_mat[idx]*this->botright->qhat[i]*scale;
            }
            this->qhat[j] = (double)(j==0?1:-1/(double)j)*this->qhat[j];
        }
#ifdef DEBUG
        if(this->level == totallevel-2){
            complex<double>* qhat = 
                (complex<double>*) malloc(p * sizeof(complex<double>));
            for(int j = 0; j<p; j++){ // row
                qhat[j] = complex<double>(0,0);
                Box* child = this->topleft;
                for(int i = 0; i < child->num; i++){// for charges in this box
                    if(j==0){// first element
                        qhat[j] += q[child->idxcharge[i]]; // how to pass q?
                    }else{ // othere elements
                        qhat[j] += ( -1./ (double) j)
                            *pow(x[child->idxcharge[i]]-this->c, j)
                            *q[child->idxcharge[i]] ;
                    }
                }
                child = this->topright;
                for(int i = 0; i < child->num; i++){// for charges in this box
                    if(j==0){// first element
                        qhat[j] += q[child->idxcharge[i]]; // how to pass q?
                    }else{ // othere elements
                        qhat[j] += ( -1./ (double) j)
                            *pow(x[child->idxcharge[i]]-this->c, j)
                            *q[child->idxcharge[i]] ;
                    }
                }
                child = this->botleft;
                for(int i = 0; i < child->num; i++){// for charges in this box
                    if(j==0){// first element
                        qhat[j] += q[child->idxcharge[i]]; // how to pass q?
                    }else{ // othere elements
                        qhat[j] += ( -1./ (double) j)
                            *pow(x[child->idxcharge[i]]-this->c, j)
                            *q[child->idxcharge[i]] ;
                    }
                }
                child = this->botright;
                for(int i = 0; i < child->num; i++){// for charges in this box
                    if(j==0){// first element
                        qhat[j] += q[child->idxcharge[i]]; // how to pass q?
                    }else{ // othere elements
                        qhat[j] += ( -1./ (double) j)
                            *pow(x[child->idxcharge[i]]-this->c, j)
                            *q[child->idxcharge[i]] ;
                    }
                }
            }
            for(int i=0; i<p; i++){
                cout << "(qtrue, qest) = ("<<qhat[i]<<","<<this->qhat[i]<<")"<<endl;
            }
            delete[] qhat;
        }
#endif
    }
}

void Box::computeincomingexp()
{
    int p = this->p;
    this->uhat = (complex<double>*) malloc(p * sizeof(complex<double>));
    // Case 1: level 0, 1
    // uhat = 0
    for(int j=0; j<p; j++){
        this->uhat[j] = complex<double>(0,0);
    }

    // Case 2: all other levels (7.4) (7.5)
    if(this->level >=2){
        // Apply T_sigma,parent^{ifi} to uhat_parent
        if(this->level >=3){
            for(int j=0; j<p;j++){//row
                for(int i=0; i<p; i++){//column
                    int idx = j*p+i;
                    this->uhat[j] += this->Tifi_mat[idx]*this->parent->uhat[i];
                }
            }
        }

        // Loop through all the boxes in the interaction list sigma
        // Apply T_sigma,interaction^{ifo} to qhat_interaction
        Node* now = this->interaction;
        while(now->data != NULL){
            for(int j=0; j<p;j++){//row
                for(int i=0; i<p; i++){//column
                    int idx = j*p+i;
                    this->uhat[j] += now->Tifo_mat[idx] * now->data->qhat[i];
                }
            }
            now = now->next;
        }
    }
}



void Box::buildTifi()
{
    //Theorem 7.2
    //c_sigma (this): (this->cx, this->cy)
    //c_tau (this->parent): (this->parent->cx, this->parent->cy)
    if(this->level < 3) return; // do nothing for box on level =0,1,2
    complex<double> cparent = this->parent->c;
    complex<double> c = this->c;

    int p = this->p;
    this->Tifi_mat = (complex<double>*) malloc(p*p * sizeof(complex<double>));
    for(int j=0; j<p; j++){ // row
        for(int i=0; i<p; i++){ //column
            int idx = j*p+i;
            if (i>=j){
                this->Tifi_mat[idx] = (double) nCr(i,j)*pow(c-cparent, i-j);
            }else{
                this->Tifi_mat[idx] = complex<double>(0,0);
            }
        }
    }
}


void Box::buildTofo()
{
    if(this->level == 0) return; // do nothing for root box
    //Theorem 7.1
    complex<double> cparent = this->parent->c;
    complex<double> c = this->c;

    int p = this->p;
    this->Tofo_mat = (complex<double>*) malloc(p*p*sizeof(complex<double>));
    for(int j=0; j<p; j++){
        for(int i=0; i<p; i++){
            int idx = j*p+i;
            if (i<=j){
                this->Tofo_mat[idx] = (double) nCr(j,i)*pow(c-cparent, j-i);
            }else{
                this->Tofo_mat[idx] = complex<double>(0,0);
            }
        }
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
    //if(this->level < 2) return; // do nothing for box on level =0,1

    Node* now = this->interaction;
    complex<double> c = this->c;
    complex<double> cinter;
    int p = this->p;

    while(now->data != NULL){
        cinter = now->data->c;
        now->Tifo_mat=(complex<double>*) malloc(p*p*sizeof(complex<double>));

        for(int j=0; j<p; j++){ // row
            for(int i=0; i<p; i++){ //column
                int idx = j*p+i; //rowise
                if (i==0) { // first column
                    if (j==0) { // first component
                        now->Tifo_mat[idx] = log(c-cinter);
                    }else{ // other components in first column
                        now->Tifo_mat[idx] = (double)(-1/(double)j)
                            *pow(cinter-c,-j); 
                    }
                }else{ // other columns
                    now->Tifo_mat[idx] = (double) pow(-1,i)* 
                        (double) nCr(i+j-1,i-1)*pow(cinter-c, -i-j);
                }
            }
        }
        now = now->next;
    }
}

void Box::buildactualpotentialbox(complex<double>* x, double* q,
        complex<double>* uapprox)
{
    // A(I_tau, I_tau)q(I_tau)
    int Nself = this->ncharge;
    int* idxself = this->idxcharge;
    for(int i=0; i<Nself; i++){
        for(int j=0; j<Nself; j++){
            if(i==j) continue;
            uapprox[idxself[i]]+=G(x[idxself[i]],x[idxself[j]])*q[idxself[j]];
        }
    }

    // A(I_tau, I_sigma)q(I_sigma) sigma is the neighbor of tau
    Box* bsigma;
    Node* now = this->neighbor;
    while(now->data != NULL){
        bsigma = now->data;
        int Nneig = bsigma->ncharge;
        int* idxneig = bsigma->idxcharge;
        for(int i=0; i<Nself; i++){
            for(int j=0; j<Nneig; j++){
                assert(idxself[i] != idxneig[j]);
                uapprox[idxself[i]]+=G(x[idxself[i]],x[idxneig[j]])*q[idxneig[j]];
            }
        }
        now = now->next;
    }

    // T^{tfi}_tau uhat_tau
#ifdef DEBUG
    int p = this->p;
    now = this->interaction;
    while(now->data != NULL){
        bsigma = now->data;
        int Nneig = bsigma->ncharge;
        int* idxneig = bsigma->idxcharge;
        for(int i=0; i<Nself; i++){
            complex<double> approx(0,0);
            complex<double> exact(0,0);

            for(int j=0; j<p; j++){
                if(j==0){
                    approx+=log(x[idxself[i]]-bsigma->c)
                        *bsigma->qhat[j];
                    uapprox[idxself[i]]+=log(x[idxself[i]]-bsigma->c)
                        *bsigma->qhat[j];
                }else{
                    approx+=pow(x[idxself[i]]-bsigma->c,-j)
                        *bsigma->qhat[j];
                    uapprox[idxself[i]]+=pow(x[idxself[i]]-bsigma->c,-j)
                        *bsigma->qhat[j];
                }
            }
            //for(int j=0; j<Nneig; j++){
            //	exact += G(x[idxself[i]],x[idxneig[j]])*q[idxneig[j]];
            //	uapprox[idxself[i]]+=G(x[idxself[i]],x[idxneig[j]])*q[idxneig[j]];
            //}
            //cout << bsigma->i << ","<<bsigma->j<<endl;
            //cout <<"Nneig"<<Nneig<<" exact = "<<exact<<", approx="<<approx<<endl;
        }
        now = now->next;
    }
#endif
    for(int i=0; i<Nself; i++){
        for(int j=0; j<this->p; j++){
            complex<double> ctau = this->c;
            //cout<<this->uhat[j]<<endl;
            uapprox[idxself[i]]+=pow(x[idxself[i]]-this->c,j) * this->uhat[j];
        }
    }
}

void Box::buildactualpotential(int totallevel, complex<double>* x, double* q,
        complex<double>* uapprox)
{
    if(this->level == totallevel-1){
        this->buildactualpotentialbox(x, q, uapprox);
        return;
    }
    this->botleft ->buildactualpotential(totallevel, x, q, uapprox);
    this->botright->buildactualpotential(totallevel, x, q, uapprox);
    this->topleft ->buildactualpotential(totallevel, x, q, uapprox);
    this->topright->buildactualpotential(totallevel, x, q, uapprox);
}
