/*----------------------------------------------------------------------------
  File    : mxFcm.c
  Contents: Compute statistics for each matrix element across multiple
            functional connectivity matrices based on input time series data.
            This file contains descriptive statistics for univariate data.
            So far, sample mean, variance, and standard deviation are
            available.
  Author  : Kristian Loewe
----------------------------------------------------------------------------*/
#include "stats.h"
#include "fcmat.h"
#include "matrix.h"
#include "edgestats.h"
#include "mex.h"

/*----------------------------------------------------------------------------
  Gateway Function
----------------------------------------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nrhs < 6 || nrhs > 7)
   ERROR("Unexpected number of input arguments.");

  // retrieve data from matlab
  REAL *data = (REAL *)mxGetData(prhs[0]);       // nodal time series data
  const mwSize *dims = mxGetDimensions(prhs[0]); // data dimensions
  DIM T = (DIM)dims[0];                          // # points in time
  DIM N = (DIM)dims[1];                          // # nodes
  int S = (int)dims[2];                          // # data sets
  int P = *(int *)mxGetData(prhs[1]);            // # threads
  int C = *(int *)mxGetData(prhs[2]);            // cache/tile size
  double maxmem = *(double *)mxGetData(prhs[3]); // max. memory
  int M = *(int *)mxGetData(prhs[4]);            // FC measure
  int f = *(int *)mxGetData(prhs[5]);            // function identifier
  DBGMSG("T: %d  N: %d  S: %d  P: %d  C: %d  maxmem: %f  M: %d  f: %d\n",
         T, N, S, P, C, maxmem, M, f);
  if (maxmem < 0)
    maxmem = 4.0;

  // create FCMAT objects
  FCMAT **fcm = malloc((size_t)S *sizeof(FCMAT*));
  for (int i = 0; i < S; i++)
    *(fcm++) = fcm_create(data+i*(T*N), N, T,
                 M|FCM_R2Z|FCM_THREAD|FCM_CACHE|FCM_MAXMEM,
                 (P != 0) ? P : 1, C, maxmem/S);
  fcm -= S;                                      // reset pointer

  // compute matrix of statistics
  if (P > 0 && C > 0 && C < N) {                 // disable multi-threading
    P = 0;                                       // if cache-based FCMAT
    DBGMSG("P: %d\n", P);                        // variant is used
  }
  MATRIX *mos = mat_create((DIM)N, (size_t)0);
  int r = 0;
  if      (f == 1)                               // 1: mean
    r = fcm_uni(fcm, S, &mean, mos, FCM_THREAD, P);
  else if (f == 2)                               // 2: std
    r = fcm_uni(fcm, S, &std,  mos, FCM_THREAD, P);
  else if (f == 3)                               // 3: var
    r = fcm_uni(fcm, S, &var,  mos, FCM_THREAD, P);
  else if (f == 4) {                             // 4: tstat2
    // grp: binary vector of length S indicating sample
    // membership:  0 -> sample #1;  1 -> sample #2
    int *grp = (int *)mxGetData(prhs[6]);
    r = fcm_tstat2(fcm, S, grp, mos, FCM_THREAD, P); }
  else if (f == 5) {                             // 5: corr
    REAL *v = (REAL *)mxGetData(prhs[6]);
    r = fcm_corr(fcm, S, v, mos, FCM_THREAD, P); }
  if (r) {                                       // clean up on error
    for (int i = 0; i < S; i++)
      fcm_delete(fcm[i]);
    free(fcm);
    mat_delete(mos);
    ERROR("Computation error.");
  }

  // copy upper triangular matrix elements to Matlab array
  mwSize E = (mwSize)N*((mwSize)N-1)/2;
  DBGMSG("E: %zu\n",E);
  if (mxIsDouble(prhs[0]))
    plhs[0] = mxCreateNumericMatrix(E, 1, mxDOUBLE_CLASS, mxREAL);
  else
    plhs[0] = mxCreateNumericMatrix(E, 1, mxSINGLE_CLASS, mxREAL);
  REAL *utmvals = (REAL *)mxGetData(plhs[0]);
  for (int i = 0; i < N; i++)
    for (int j = i+1; j < N; j++)
      *utmvals++ = mat_get(mos, (DIM)i, (DIM)j);

  // clean up
  for (int i = 0; i < S; i++)
    fcm_delete(fcm[i]);
  free(fcm);
  mat_delete(mos);

}  // mexFunction()
