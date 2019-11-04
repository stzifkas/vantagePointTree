#include "vptree.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double *distance_array;
int array_index;

void  *generateInner(void *arg);
void * generateOuter(void *arg);
vptree * generatevp(vptree *T);
double calculate_md(vptree * T);


double qselect(double *v, int len, int k)
{
#	define SWAP(a, b) { tmp = v[a]; v[a] = v[b]; v[b] = tmp; }
  int i, st;
  double tmp;
	for (st = i = 0; i < len - 1; i++) {
		if (v[i] > v[len-1]) continue;
		SWAP(i, st);
		st++;
	}

	SWAP(len-1, st);

	return k == st	?v[st]
			:st > k	? qselect(v, st, k)
				: qselect(v + st, len - st, k - st);
}

vptree * generatevp(vptree *T)
{
  if ((T->n) == 0)
  {
    T = NULL;
    return T;
  }
  else if ((T->n) == 1)
  {
    return T;
  }
  int starting_index = array_index;
  T->starting_index = array_index; //T will start at starting_index
  double * vp = getVP(T); //get VP;
  int i;
  for (i=0;i<((T->n)-1);i++)
  {
    double temp = 0;
    int j;
    for(j=0;j<T->d;j++)
    {
      temp = temp + ((*(T->points + i*T->d + j)) - (*(vp + j)))*((*(T->points + i*T->d + j)) - (*(vp + j)));
    }
    distance_array[array_index] = sqrt(temp);
    T->ending_index = array_index;
    array_index++;
  }
  T->md = calculate_md(T);
  generateInner((void *)T);
  generateOuter((void *)T);
  vptree *ret_inner;
  vptree *ret_outer;
  T->outer = generatevp(T->outer);
  T->inner = generatevp(T->inner);
  return T;
}

vptree * getInner(vptree * T)
{
  return T->inner;
}

vptree * getOuter(vptree * T)
{
  return T->outer;
}

vptree * buildvp(double *X, int n, int d)
{
  array_index = 0;
  int size = ((log(n)/log(2))*n); //log(n)/log(2) is tree height, n is distances width
  distance_array = calloc(size,sizeof(double));
  int *idx = malloc(n*sizeof(int)); //will hold initial indexes
  int i;
  for (i=0;i<n;i++)
  {
    idx[i] = i; //populate with initial indexes
  }
  vptree *T = malloc(sizeof(vptree));
  T->points = X;
  T->n = n;
  T->d = d;
  T->idxs = idx;
  T->inner = NULL;
  T->outer = NULL;
  T = (vptree *) generatevp((void *)T);
  return T;
}
void * generateInner(void * arg)
{
  vptree * T = (vptree *) arg;
  double median_distance = getMD(T);
  int sz = (T->n);
  vptree * inner_tree = malloc(sizeof(vptree));
  double * inner = malloc(sz*(T->d)*sizeof(double)); //allocate on maximum (impossible) size, fix it later
  int * idx = malloc(sz*sizeof(int)); //allocate on maximum (impossible) size, fix it later
  int inner_n = 0;
  int i,j;
  int n = T->n;
  int d = T->d;
  int *idxs = T->idxs;
  double *vp = getVP(T);
  double *points = T->points;
  for (i=T->starting_index;i<=T->ending_index;i++)
  {
    if (distance_array[i]<=median_distance)
    {
      inner_n++;
      idx[inner_n-1] = T->idxs[i-T->starting_index];
      for (j=0;j<d;j++)
      {
        *(inner + ((inner_n-1)*d+j)) = *(points + ((i-T->starting_index)*d + j));
      }
    }
  }
  inner_tree->n = inner_n;
  inner = realloc(inner,(inner_n)*d*sizeof(double));
  idx = realloc(idx,(inner_n)*sizeof(int));
  inner_tree->d = d;
  inner_tree->idxs = idx;
  inner_tree->points = inner;
  inner_tree->inner = NULL;
  inner_tree->outer = NULL;
  inner_tree->md = 0;
  T->inner = inner_tree;
  return NULL;
}

void * generateOuter(void *arg)
{
  vptree * T = (vptree *) arg;
  double median_distance = getMD(T);
  vptree * outer_tree = malloc(sizeof(vptree));
  int i,j;
  int n = T->n;
  int d = T->d;
  int sz = (T->n);
  double *outer = malloc(sz*d*sizeof(double)); //will hold outer subtree elements
  int *idx = malloc(sz*sizeof(int)); //will hold initial indexes of outer elements
  int outer_n = 0;
  double *vp = getVP(T);
  double *points = T->points;
  for (i=T->starting_index;i<=T->ending_index;i++)
  {
    if (distance_array[i]>median_distance)
    {
      outer_n++;
      idx[outer_n-1] = T->idxs[i-T->starting_index]; //populate idx
      for (j=0;j<d;j++)
      {
        *(outer + ((outer_n-1)*d+j)) = *(points + ((i-T->starting_index)*d + j));
      }
    }
  }
  outer = realloc(outer,(outer_n)*d*sizeof(double)); //fix size
  idx = realloc(idx,(outer_n)*sizeof(int)); //fix size
  outer_tree->n = outer_n; //outer_n -> number of outer elements
  outer_tree->idxs = idx;
  outer_tree->points = outer;
  outer_tree->d = d;
  outer_tree->inner = NULL;
  outer_tree->outer = NULL;
  outer_tree->md = 0;
  T->outer = outer_tree;
}

double calculate_md(vptree * T)
{
  int k;
  int length = T->n - 1;
  if (length==0)
    return 0;
  if (length%2 != 0)
    k = (length)/2 ; //median is Kth element
  else
    k = (length/2)-1; //median is Kth element
  double distanceCopy[length]; //qselect modifies array so we have to copy it
  int i;
  int j=0;
  for (i=T->starting_index;i<=T->ending_index;i++)
  {
    distanceCopy[j] = distance_array[i];
    j++;
  }
  return qselect(distanceCopy,length,k);
}

double getMD(vptree *T)
{
  return T->md;
}

double * getVP(vptree *T)
{
  double *vp = malloc((T->d)*sizeof(double));
  int i;
  for (i=0;i<(T->d);i++)
  {
    *(vp+i) = *(T->points + (T->d)*((T->n)-1)+i);
  }
  return vp;
}

int getIDX(vptree *T)
{
  int idx;
  if (T->n>0)
    idx = T->idxs[(T->n)-1];
  else
    idx = -1;
  return idx;
  }
