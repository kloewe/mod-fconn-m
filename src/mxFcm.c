/*----------------------------------------------------------------------------
  File    : mxFcm.c
  Contents: Compute statistics for each matrix element across multiple
            functional connectivity matrices based on input time series data.
            This file contains descriptive statistics for univariate data.
            So far, sample mean, variance, and standard deviation are
            available.
  Author  : Kristian Loewe
----------------------------------------------------------------------------*/
#include <string.h>
#include "stats.h"
#include "fcmat.h"
#include "matrix.h"
#include "edgestats.h"
#include "mex.h"

#ifndef REAL_IS_DOUBLE
#  error "REAL_IS_DOUBLE is undefined."
#else
#  if REAL_IS_DOUBLE
#    define REAL_CN double                         // matlab class name
#  else
#    define REAL_CN single
#  endif
#endif

/*----------------------------------------------------------------------------
  Array of Func structs
----------------------------------------------------------------------------*/
static Func functions[12] = {
  { .kind = 1, .u.f1 = &sum_w },            // sum      1
  { .kind = 1, .u.f1 = &mean_w },           // mean     2
  { .kind = 1, .u.f1 = &var_w },            // var      3
  { .kind = 1, .u.f1 = &std_w },            // std      4
  { .kind = 1, .u.f1 = &tstat_w },          // tstat    5
  { .kind = 1, .u.f1 = &mdiff_w },          // mdiff    6
  { .kind = 1, .u.f1 = &tstat2_w },         // tstat2   7
  { .kind = 1, .u.f1 = &pairedt_w },        // pairedt  8
  { .kind = 1, .u.f1 = &didt_w },           // didt     9
  { .kind = 2, .u.f2 = &corrv_w },          // corrv   10
};

/*----------------------------------------------------------------------------
  Gateway Function
----------------------------------------------------------------------------*/

// mxFcm(f, data, n12, ..., p, c, mm, fcmeas) 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nrhs < 3 || nrhs > 10)
    ERROR("Unexpected number of input arguments: %d", nrhs);

  // retrieve data from matlab
  if (!mxIsClass(prhs[0], "int32"))                // function identifier
    ERROR("Argument %d must be of type '%s'.", 0, "int32");
  int fno = *(int *)mxGetData(prhs[0]);

  if (!mxIsClass(prhs[1], TOSTRING(REAL_CN)))      // data
    ERROR("Argument %d must be of type '%s' but is of type '%s'", 1,
            TOSTRING(REAL_CN), mxGetClassName(prhs[1]));
  REAL *data = (REAL *)mxGetData(prhs[1]);

  const mwSize *dims = mxGetDimensions(prhs[1]);   // data dimensions
  int T = (int)dims[0];                            // # points in time
  int N = (int)dims[1];                            // # nodes
  int S = (int)dims[2];                            // # data sets

  int *AS = (int *)mxGetData(prhs[2]);             // [S1 S2 ... ]

  if (!mxIsClass(prhs[nrhs-4], "int32"))           // # threads
    ERROR("Argument %d must be of type '%s' but is of type '%s'.", nrhs-4,
          "int32", mxGetClassName(prhs[nrhs-4]));
  int P = *(int *)mxGetData(prhs[nrhs-4]);

  if (!mxIsClass(prhs[nrhs-3], "int32"))           // cache/tile size
    ERROR("Argument %d must be of type '%s' but is of type '%s'.", nrhs-3,
          "int32", mxGetClassName(prhs[nrhs-3]));
  int C = *(int *)mxGetData(prhs[nrhs-3]);
  
  if (!mxIsClass(prhs[nrhs-2], TOSTRING(REAL_CN))) // max. memory
    ERROR("Argument %d must be of type '%s' but is of type '%s'.", nrhs-2,
          TOSTRING(REAL_CN), mxGetClassName(prhs[nrhs-2]));
  REAL mm = *(REAL *)mxGetData(prhs[nrhs-2]);
  
  if (!mxIsClass(prhs[nrhs-1], "int32"))           // FC measure
    ERROR("Argument %d must be of type '%s' but is of type '%s'.", nrhs-1,
            "int32", mxGetClassName(prhs[nrhs-1]));
  int M = *(int *)mxGetData(prhs[nrhs-1]);

  DBGMSG("T: %d  N: %d  S: %d  P: %d  C: %d  mm: %f  M: %d  fno: %d\n",
         T, N, S, P, C, mm, M, fno);

  if (mm < 0)                                      // max. memory
    mm = 4.0;                                      // (default choice)

  DBGMSG("T: %d  N: %d  S: %d  P: %d  C: %d  mm: %f  M: %d  fno: %d\n",
          T, N, S, P, C, mm, M, fno);

  REAL *v = NULL;                                  // optional add. argument
  REAL *a = NULL;                                  // pre-normalized add. arg.
  void *mem = NULL;                                // mem. for pre-norm. a.a.

  if (fno == 10) {                                 // function is 'corrv'
    v = (REAL *)mxGetData(prhs[3]);

    // get aligned memory for the pre-normalized values
    mem = malloc((size_t)S *sizeof(REAL) +31);
    if (!mem)
      ERROR("ERROR: malloc failed");
    a = (REAL*)(((uintptr_t)mem +31) & ~(uintptr_t)31);

    // pre-normalize the add. variable
    REAL sqr = 0;
    REAL ma = mean(v, S);
    for (int k = 0; k < S; k++) {
      a[k] = v[k] - ma;
      sqr += a[k]*a[k]; }
    sqr = sqrt(sqr);
    sqr = (sqr > 0) ? 1/sqr : 0;
    for (int k = 0; k < S; k++)
      a[k] *= sqr;
  }

  // create FCMAT objects
  FCMAT **fcm = malloc((size_t)S *sizeof(FCMAT*));
  if (!fcm) ERROR("Malloc failed.");

  for (int i = 0; i < S; i++) {
    fcm[i] = fcm_create(data+i*(T*N), N, T,
                          M|FCM_R2Z|FCM_THREAD|FCM_CACHE|FCM_MAXMEM,
                          (P != 0) ? P : 1, C, mm/(REAL)S);
    if (!fcm[i]) ERROR("Could not create FCMAT. i: %d", i); }

  // disable multi-threading if cache-based FCMAT variant is used
  if (P > 0 && C > 0 && C < N) {
    P = 0;
    DBGMSG("P: %d\n", P);
  }

  // init matrix of statistics
  MATRIX *mos = mat_create((DIM)N, (size_t)0);
  if (!mos) ERROR("Could not create MATRIX.");

  // compute statistics
  int r = fcm_uni(fcm, S, functions[fno-1], AS, a, mos, FCM_THREAD, P);

  // free mem if necessary
  if (fno == 10)
    free(mem);

  if (r) {                                       // clean up on error
    for (int i = 0; i < S; i++)
      fcm_delete(fcm[i]);
    free(fcm);
    mat_delete(mos);
    ERROR("Computation error.");
  }

  // copy upper triangular matrix elements to Matlab array
  mwSize E = (mwSize)N*((mwSize)N-1)/2;
  DBGMSG("%s\n", TOSTRING(REAL));
  DBGMSG("fno: %d\n", fno);
  DBGMSG("E: %zu\n",E);
  #if REAL_IS_DOUBLE
  plhs[0] = mxCreateNumericMatrix(E, 1, mxDOUBLE_CLASS, mxREAL);
  #else
  plhs[0] = mxCreateNumericMatrix(E, 1, mxSINGLE_CLASS, mxREAL);
  #endif
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
