/* GS smoothing to solve -Laplace u =f in 2D
 * Author: Xinyang Wang
 */
#include <stdio.h>
#include <math.h>
#include "util.h"
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif


/* compuate global residual, assuming ghost values are updated */
double compute_residual(double **u, int N, double invhsq)
{
  int i, j;
  double tmp, res = 0.0;
  #pragma omp parallel for default(none) shared(u,N,invhsq) private(i,j,tmp) reduction(+:res)
  for (i = 1; i <= N; i++) {
  	for (j = 1; j <= N; j++) {
    	tmp = ((4.0*u[i][j] - u[i-1][j] - u[i][j-1]- u[i+1][j] - u[i][j+1]) * invhsq - 1);
    	res += tmp * tmp;
  	}
  }
  
  return sqrt(res);
}


int main(int argc, char * argv[])
{
  int i, j, N, iter, max_iters;

  sscanf(argv[1], "%d", &N);
  sscanf(argv[2], "%d", &max_iters);

#ifdef _OPENMP
#pragma omp parallel
  {
    int my_threadNum = omp_get_thread_num();
    int numThreads = omp_get_num_threads();
    printf("Hello from thread %d of %d\n",my_threadNum, numThreads);
  }
#else
  printf("Hello from process %d of %d\n",0, 1);
#endif

  
  /* timing */
  timestamp_type time1, time2;
  get_timestamp(&time1);

  /* Allocation of matrix, including boundary ghost points */
  double ** u    = (double **) calloc(sizeof(double *), N+2);
  	for (i = 0; i < N+2; i++) {
  		u[i] = (double *) calloc(sizeof(double), N+2);
  	}

  double h = 1.0 / (N + 1);
  double hsq = h * h;
  double invhsq = 1./hsq;
  double res, res0, tol = 1e-5;

  /* initial residual */
  res0 = compute_residual(u, N, invhsq);
  res = res0;
  for (i=0; i < N+2 ; i++) {
	u[i][0] = u[i][N+1] = 0.0;
	u[0][i] = u[N+1][i] = 0.0;
  }

  for (iter = 0; iter < max_iters && res/res0 > tol; iter++) {

# pragma omp parallel for default(none) shared (hsq,N,u) private(i,j)
    /* update black points */
  	for (i = 1; i <= N; i++) {     
  	  for (j = (i%2)+1; j <= N; j+=2) {     
    //	for (j = 1; j <= N; j++) {
    //		if ( ( i + j ) % 2 == 0) {
      			u[i][j]  = 0.25 * (hsq + u[i-1][j] + u[i][j-1] + u[i+1][j] + u[i][j+1]);
     // 		}
      }
    }
# pragma omp parallel for default(none) shared (hsq,N,u) private(i,j)
    /* update red points */
  	for (i = 1; i <= N; i++) {
      for (j = (i+1)%2+1; j <= N; j+=2) {  	
   // 	for (j = 1; j <= N; j++){
   // 	    if ( ( i + j ) % 2 == 1) {
      		u[i][j]  = 0.25 * (hsq + u[i-1][j] + u[i][j-1] + u[i+1][j] + u[i][j+1]);
    //  		}
      }
    }
    if (0 == (iter % 10)) {
      res = compute_residual(u, N, invhsq);
      printf("Iter %d: Residual: %g\n", iter, res);
    }
  }

  /* Clean up */
  free(u);

  /* timing */
  get_timestamp(&time2);
  double elapsed = timestamp_diff_in_seconds(time1,time2);
  printf("Time elapsed is %f seconds.\n", elapsed);
  return 0;
}
