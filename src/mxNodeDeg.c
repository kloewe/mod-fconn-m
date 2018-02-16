/*----------------------------------------------------------------------------
  File    : mxNodeDeg.c
  Contents: Compute the degree of each node in a functional connectivity graph
  Author  : Kristian Loewe
----------------------------------------------------------------------------*/
#include "fcmat.h"
#include "nodedeg.h"
#include "mex.h"

/*----------------------------------------------------------------------------
  Gateway Function
----------------------------------------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nrhs != 6)
    ERROR("Unexpected number of input arguments.");

  // retrieve data from matlab
  REAL *data = (REAL *)mxGetData(prhs[0]);       // nodal time series data
  const mwSize *dims = mxGetDimensions(prhs[0]); // data dimensions
  DIM T = (DIM)dims[0];                          // # points in time
  DIM N = (DIM)dims[1];                          // # nodes
  REAL thr = *(REAL *)mxGetData(prhs[1]);        // FC threshold
  int P = *(int *)mxGetData(prhs[2]);            // # threads
  int C = *(int *)mxGetData(prhs[3]);            // cache/tile size
  double maxmem = *(double *)mxGetData(prhs[4]); // max. memory
  int M = *(int *)mxGetData(prhs[5]);            // FC measure
  DBGMSG("T: %d  N: %d  P: %d  C: %d  maxmem: %f  M: %d  thr: %f\n",
          T, N, P, C, maxmem, M, thr);

  // memory allocation for result
  plhs[0] = mxCreateNumericMatrix((mwSize)N, 1, mxINT32_CLASS, mxREAL);
  int *intres = (int *)mxGetData(plhs[0]);

  // create FCMAT object
  FCMAT *fcm = fcm_create(data, N, T, M|FCM_THREAD|FCM_CACHE|FCM_MAXMEM,
                          (P != 0) ? P : 1, C, maxmem);
  if (!fcm)
    ERROR("Memory allocation failed.");

  // degree computation
  if (P > 0 && C > 0 && C < N) {                 // disable multi-threading
    P = 0;                                       // if cache-based FCMAT
    DBGMSG("P: %d\n", P);                        // variant is used
  }
  int r = fcm_nodedeg(fcm, thr, intres, FCM_THREAD, P);
  if (r) {                                       // clean up on error
    fcm_delete(fcm);
    ERROR("Computation error.");
  }

  // clean up
  fcm_delete(fcm);

}  // mexFunction()
