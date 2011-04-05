
#include <mex.h>
#include "_reg_array_interface.h"

#include <limits>
#include <string.h>
#include <math.h>
#include <cmath>

/*###############################################################################################*/
/* Matlab extensions */
/*###############################################################################################*/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   int status = 1;

   /* Check for correct number of arguments. */
   if (nrhs!=3){
       mexErrMsgTxt("3 inputs required: Affine transformation [4,4], image size [1,3], control points grid pacing [1,3].");
   }
   if (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS || mxGetClassID(prhs[1]) != mxDOUBLE_CLASS || mxGetClassID(prhs[2]) != mxDOUBLE_CLASS)
       mexErrMsgTxt("Inputs must be real doubles.");

   const mxClassID cid_affine      = mxGetClassID(prhs[0]);
   const int       dim_affine      = mxGetNumberOfDimensions(prhs[0]); 

   const mxClassID cid_image_size  = mxGetClassID(prhs[1]);
   const int       dim_image_size  = mxGetNumberOfDimensions(prhs[1]); 

   const mxClassID cid_cp_spacing  = mxGetClassID(prhs[2]);
   const int       dim_cp_spacing  = mxGetNumberOfDimensions(prhs[2]); 

   //extract pointers to input matrices
   double* affineTransformation = (double *) (mxGetData(prhs[0]));
   double* imageSize            = (double *) (mxGetData(prhs[1]));
   double* gridSpacing          = (double *) (mxGetData(prhs[2]));

   //convert data types
   float affineTransformation_f[16]; for (int i=0; i<16; i++) affineTransformation_f[i] = affineTransformation[i];
   int imageSize_i[3]; imageSize_i[0] = imageSize[0]; imageSize_i[1] = imageSize[1]; imageSize_i[2] = imageSize[2];
   float gridSpacing_f[3]; gridSpacing_f[0] = gridSpacing[0]; gridSpacing_f[1] = gridSpacing[1]; gridSpacing_f[2] = gridSpacing[2];

   //calculate size of control points grid, in order to allocate output
   int cp_size[3];
   reg_array_compute_control_point_size(cp_size, imageSize_i, gridSpacing_f);

   //allocate output matrix for control points   
   mwSize mw_cp_size[4];
   mw_cp_size[0] = (mwSize)cp_size[0];
   mw_cp_size[1] = (mwSize)cp_size[1];
   mw_cp_size[2] = (mwSize)cp_size[2];
   mw_cp_size[3] = 4;

   plhs[0] =  mxCreateNumericArray(4, mw_cp_size, mxSINGLE_CLASS, mxREAL);
   float *cp_ptr = (float *)(mxGetData(plhs[0]));

   /* Compute control points from affine */
   status = reg_array_bspline_initialiseControlPointGridWithAffine(cp_ptr, affineTransformation_f, imageSize_i, gridSpacing_f);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while performing Affine Transform.");
   return;
}


