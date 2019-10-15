#include "mex.h"
#include "mexBistocNormalize_match.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    enum{X, idx1, ID1, idx2, ID2, tol};
    enum{Y};
       
    int N;
    
    double* pX = mxGetPr(prhs[X]);
    int* pIdx1 = (int*)mxGetPr(prhs[idx1]);
    int* pID1 = (int*)mxGetPr(prhs[ID1]);
    int* pIdx2 = (int*)mxGetPr(prhs[idx2]);
    int* pID2 = (int*)mxGetPr(prhs[ID2]);
    double* pTol = mxGetPr(prhs[tol]);
    
    N = mxGetM(prhs[X]);
        
    plhs[Y] = mxCreateDoubleMatrix(N, 1, mxREAL);
    double* pY = mxGetPr(plhs[Y]);
            
    bistocNorm(pX, N, pIdx1, pID1, pIdx2, pID2, pTol, pY);
}