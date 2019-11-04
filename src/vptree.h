#ifndef VPTREE_H
#define VPTREE_H

typedef struct vptree{
    double *points; //current dataset
    int *idxs; //initial indexes for current points
    int n; //number of current points
    int d; //points dimensions
    int starting_index; //here starts distance_array part of T
    int ending_index; //last element of distance_array part of T
    struct vptree *inner; //pointer to inner subtree
    struct vptree *outer; //pointer to outer subtree
    double md; //median distance
}vptree;

vptree * buildvp(double *X, int n, int d);
vptree * getInner(vptree *T);
vptree * getOuter(vptree *T);
double getMD(vptree *T);
double * getVP(vptree *T);
int getIDX(vptree *T);
#endif
