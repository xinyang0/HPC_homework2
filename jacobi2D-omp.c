/* Jacobi smoothing to solve -Laplace u=f in 2D
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
# pragma omp parallel
  {
#ifdef _OPENMP
    int my_threadnum = omp_get_thread_num();
    int numthreads = omp_get_num_threads();
#else
    int my_threadnum = 0;
    int numthreads = 1;
#endif
    printf("Hello, I'm thread %d out of %d\n", my_threadnum, numthreads);
  }

  /* timing */
  timestamp_type time1, time2;
  get_timestamp(&time1);

  /* Allocation of matrix, including boundary ghost points */
  double ** u    = (double **) calloc(sizeof(double *), N+2);
  	for (i = 0; i < N+2; i++) {
  		u[i] = (double *) calloc(sizeof(double), N+2);
  	}
  
  double ** unew    = (double **) calloc(sizeof(double *), N+2);
  	for (i = 0; i < N+2; i++) {
  		unew[i] = (double *) calloc(sizeof(double), N+2);
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

#pragma omp parallel for default(none) shared(N,unew,u,hsq) private(i,j)
    /* Jacobi step for all the inner points */
    for (i = 1; i <= N; i++) {
    	for (j = 1; j <= N; j++){
      		unew[i][j]  = 0.25 * (hsq + u[i-1][j] + u[i][j-1] + u[i+1][j] + u[i][j+1]);
      }
    }

    /* copy new_u onto u */
    double **utemp;
    utemp = u;
    u = unew;
    unew = utemp;
    //    memcpy(u, unew, (N+2)*sizeof(double));
    if (0 == (iter % 100)) {
      res = compute_residual(u, N, invhsq);
      printf("Iter %d: Residual: %g\n", iter, res);
    }
  }

  /* Clean up */
  	for (i = 0; i < N+2; i++) {
  		free(u[i]);
  	}
  	free(u);
    
    for (i = 0; i < N+2; i++) {
  		free(unew[i]);
  	}
  	free(unew);

  /* timing */
  get_timestamp(&time2);
  double elapsed = timestamp_diff_in_seconds(time1,time2);
  printf("Time elapsed is %f seconds.\n", elapsed);
  return 0;
}
