#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>


int main(void){
    int i, j, k;
    struct timespec start, stop; 
    double time;
    //int n = 4096; // matrix size is n*n
    int n = 4096;
    double **A = (double**) malloc (sizeof(double*)*n);
    double **B = (double**) malloc (sizeof(double*)*n);
    double **C = (double**) malloc (sizeof(double*)*n);
    for (i=0; i<n; i++) {
      A[i] = (double*) malloc(sizeof(double)*n);
      B[i] = (double*) malloc(sizeof(double)*n);
      C[i] = (double*) malloc(sizeof(double)*n);
    }
    
    for (i=0; i<n; i++){
      for(j=0; j< n; j++){
        A[i][j]=i;
        B[i][j]=i+j;
        C[i][j]=0;      
      }
    }
        
    if( clock_gettime(CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}
    
    // Your code goes here //
    // Matrix C = Matrix A * Matrix B //  
    //*******************************//
    #pragma omp parallel for
    for(i=0; i< n; i++) {
      for (j=0; j< n; j++) {
        for (k=0; k< n; k++) {
          C[i][j] = C[i][j] + A[i][k] * B[k][j]; 
        }
      }
    }
    
    //*******************************//
    
    if( clock_gettime( CLOCK_REALTIME, &stop) == -1 ) { perror("clock gettime");}   
    time = (stop.tv_sec - start.tv_sec)+ (double)(stop.tv_nsec - start.tv_nsec)/1e9;
    
    printf("Number of FLOPs = %lu, Execution time = %f sec,\n%lf MFLOPs per sec\n", 2*n*n*n, time, 1/time/1e6*2*n*n*n);   
    printf("C[100][100]=%f\n", C[100][100]);
    
    // release memory
    for (i=0; i<n; i++) {
      free(A[i]);
      free(B[i]);
      free(C[i]);
    }
    free(A);
    free(B);
    free(C);
    return 0;
}



#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>

int main(void){          // this example is about multiplying a 4096*4096 matrix with a vector
  int i, j;
  int n = 4096;
  struct timespec start, stop;
  double time;

  double **matrix = (double**) malloc (sizeof(double*)*n);
  for (i=0; i<n; i++) matrix[i] = (double*) malloc(sizeof(double)*n);

  double * vector = (double *) malloc (sizeof(double)*n);
  double * result = (double *) malloc (sizeof(double)*n);
  for (i=0; i<n; i++){
  	for (j=0; j<n; j++) matrix[i][j]=2;
  	vector[i]=3;
  	result[i]=0;
  }
  if( clock_gettime( CLOCK_REALTIME, &start) == -1 ) { perror( "clock gettime" );}

  // omp_set_num_threads(8);

  // matrix * vector
  #pragma omp parallel for
  for(i=0; i< n; i++)
    for (j=0; j< n; j++)
      result[i] = result[i] + matrix[i][j] * vector[j];


  if( clock_gettime( CLOCK_REALTIME, &stop) == -1 ) { perror( "clock gettime" );}
  time = (stop.tv_sec - start.tv_sec)+ (double)(stop.tv_nsec - start.tv_nsec)/1e9;

  printf("Number of FLOPs = %d, Execution time = %f sec,\n%f MFLOPs per sec\n", 2*n*n, time, 2*n*n/time/1e6);

  // release memory
  free(vector);
  free(result);
  for (i=0; i<n; i++)
  	free(matrix[i]);

  return 0;
}