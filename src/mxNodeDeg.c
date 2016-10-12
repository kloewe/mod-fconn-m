/*----------------------------------------------------------------------------
  File    : mxNodedeg.c
  Contents: Compute the degree of each node in a functional connectivity graph
  Author  : Kristian Loewe
----------------------------------------------------------------------------*/
#include "cpuinfo.h"
#include "fcmat.h"
#include "nodedeg.h"
#include "mex.h"

/*----------------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------------*/
#ifndef NDEBUG
#  ifndef __STRING
#    define __STRING(x) #x
#  endif
#  define TOSTRING(x) __STRING(x)
#  define ERROR(s) mexErrMsgTxt(__FILE__ ":" TOSTRING(__LINE__) ": "s)
#else
#  define ERROR(s) mexErrMsgTxt((s))
#endif

/*----------------------------------------------------------------------------
  Gateway Function
----------------------------------------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nrhs != 5)
    ERROR("Unexpected number of input arguments.");
  
  const mwSize *dims = mxGetDimensions(prhs[0]); // dimensions of input data
  REAL thr = *(REAL *)mxGetData(prhs[1]);        // FC threshold
  DIM T = (DIM)dims[0];                          // # points in time
  DIM N = (DIM)dims[1];                          // # nodes
  REAL *data = (REAL *)mxGetData(prhs[0]);       // nodal time series
  FCMAT *fcm;                                    // FC matrix
  int mode = *(int *)mxGetData(prhs[2]);         // FC measure
  #ifndef NDEBUG
  printf("T: %d  N: %d  thr: %f", T, N, thr); mexEvalString("drawnow");
  #endif

  // memory allocation for result
  plhs[0] = mxCreateNumericMatrix((mwSize)N, 1, mxINT32_CLASS, mxREAL);
  int *intres = (int *)mxGetData(plhs[0]);
  
  // tile size
  int C = (size_t) *(int *)mxGetData(prhs[4]);
  if (C == -1) {   // auto-determine
    #ifndef NDEBUG
    printf("  C: %d", C); mexEvalString("drawnow");
    #endif
    if (N < 20000)
      C = N;       // half-stored
    else {
      if (corecnt() > 6)
        C = 1024;  // cache-based (tile size: 1024)
      else
        C = 0;     // on-demand
    }
    #ifndef NDEBUG
    printf(" => %d", C); mexEvalString("drawnow");
    #endif
  } else {
    #ifndef NDEBUG
    printf("  C: %d", C); mexEvalString("drawnow");
    #endif
  }
  if ((C < 0) || (C > N))
    ERROR("Invalid tile size specified.");
  
  // determine number of threads and start degree computation
  int rc;
  int P = *(int *)mxGetData(prhs[3]);
  #ifndef NDEBUG
  printf("  P: %d\n", P); mexEvalString("drawnow");
  #endif
  if (P == -1) { // auto-determine
    // FCMAT's number of threads will be auto-determined by FCMAT
    fcm = fcm_create(data, N, T, mode|FCM_CACHE, C);
    if (!fcm) ERROR("Memory allocation failed.");
    // the number of threads to be used for degree computation will be
    // auto-determined by fcm_nodedeg for C == 0 (on-demand) and C == N
    // half-stored) but for 0 < C < N (cache-based) no multi-threading
    // can be used by fcm_nodedeg.
    if (C == 0 || C == N)
      rc = fcm_nodedeg(fcm, thr, intres, 0);
    else
      rc = fcm_nodedeg(fcm, thr, intres, FCM_THREAD, 0); }
  else {         // use specified value
    // FCMAT's number of threads is set to the specified value
    fcm = fcm_create(data, N, T, mode|FCM_THREAD|FCM_CACHE, C, P);
    // the number of threads is set to the specified value unless the
    // cache-based FCMAT variant is to be used in which case no
    // multi-threading can be used by fcm_nodedeg
    if (C == 0 || C == N)
      rc = fcm_nodedeg(fcm, thr, intres, FCM_THREAD, P);
    else
      rc = fcm_nodedeg(fcm, thr, intres, FCM_THREAD, 0);
  }
  if (rc)
    ERROR("Computation error.");
  
  fcm_delete(fcm);
}  // mexFunction()
