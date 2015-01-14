#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdarg.h>
#define FLOAT double
 
static int compar (const void *a, const void *b);
void sort(int n, FLOAT *arr, int *idx);
void permutation_FLOAT(int array_size, int *idx, int count, ...);
void permutation_int(int array_size, int *idx, int count, ...);


FLOAT *base_arr;
static int compar (const void *a, const void *b)
{
  int aa = *((int *) a), bb = *((int *) b);
  if (base_arr[aa] < base_arr[bb])
    return -1;
  if (base_arr[aa] == base_arr[bb])
    return 0;
  if (base_arr[aa] > base_arr[bb])
    return 1;
}
 
void sort(int n, FLOAT *arr, int *idx)
{
  int i,j;
  FLOAT *copy = malloc (sizeof (FLOAT) * n);
  for (i = 0; i < n; i++)
    {
      idx[i] = i;
    }
  base_arr = arr;
  qsort (idx, n, sizeof (int), compar);
  for(j=0;j<n;j++)
  {
      copy[j]=arr[j];
  }
  for(j=0;j<n;j++)
  {
      arr[j]=copy[idx[j]];
  }
}


void permutation_FLOAT(int array_size, int *idx, int count, ...)
{
    va_list ap;
    int i,j;
    va_start (ap, count);         /* Initialize the argument list. */
    for (i = 0; i < count; i++)
    {
        FLOAT *array = va_arg (ap, FLOAT*);    /* Get the next argument value. */
        FLOAT *copy = malloc(array_size * sizeof(FLOAT));
        for(j=0;j<array_size;j++)
            copy[j]=array[j];
        for(j=0;j<array_size;j++)
            array[j]=copy[idx[j]];
    }   
    va_end (ap);                  /* Clean up. */
}

void permutation_int(int array_size, int *idx, int count, ...)
{
    va_list ap;
    int i,j;
    va_start (ap, count);         /* Initialize the argument list. */
    for (i = 0; i < count; i++)
    {
        int *array = va_arg (ap, int*);    /* Get the next argument value. */
        int *copy = malloc(array_size * sizeof(int));
        for(j=0;j<array_size;j++)
            copy[j]=array[j];
        for(j=0;j<array_size;j++)
            array[j]=copy[idx[j]];
    }   
    va_end (ap);                  /* Clean up. */
}



int main()
{

    /* int i,N=10; */
    /* FLOAT *x = malloc(N * sizeof(FLOAT)); */
    /* FLOAT *y = malloc(N * sizeof(FLOAT)); */
    /* FLOAT *z = malloc(N * sizeof(FLOAT)); */
    /* int *idx = malloc(N * sizeof(int)); */
    
    /* srand(time(NULL)); */
    /* for (i=0;i<N;i++) */
    /* {   x[i] = rand(); */
    /*     y[i] = i; */
    /*     z[i]=i*i; */
    /* } */
    
    /* for (i=0;i<N;i++) */
    /* { */
    /*     printf("%f %f  %f\n",x[i],y[i],z[i]); */
    /* } */

    /* sort(N,x,idx); */
    /* permutation(N,idx,2,y,z); */
    /* printf("Sorted.\n"); */

    /* for (i=0;i<N;i++) */
    /* { */
    /*     printf("%f %f  %f\n",x[i],y[i],z[i]); */
    /* } */


}



