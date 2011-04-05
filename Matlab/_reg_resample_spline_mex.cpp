


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
   /* Check for proper number of arguments. */
   if (nrhs!=3){
      mexErrMsgTxt("3 inputs required: Image [N,M,K], Spline Nodes [Y,W,Z,4], Nodes Spacing [y,w,z]");
   }

   int   image_size[3];     //
   float spacing[3];        //
   int   nodes_size[3];
   int   enable_gpu;        // Flag for GPU Acceleration: 1 to enable, 0 to disable.
   int   status = 1;        // Return 0 if everython ok, 1 if errors.

   image_size[0] = mxGetDimensions(prhs[0])[0];
   image_size[1] = mxGetDimensions(prhs[0])[1];
   image_size[2] = mxGetDimensions(prhs[0])[2];

   nodes_size[0] = mxGetDimensions(prhs[1])[0];
   nodes_size[1] = mxGetDimensions(prhs[1])[1];
   nodes_size[2] = mxGetDimensions(prhs[1])[2];

   spacing[0] = ((float*) mxGetData(prhs[2]))[0];
   spacing[1] = ((float*) mxGetData(prhs[2]))[1];
   spacing[2] = ((float*) mxGetData(prhs[2]))[2];

   /* The inputs must be noncomplex floating point matrices */
   for (int n=0; n<2; n++)
      {
      if ( mxGetClassID(prhs[n]) != mxSINGLE_CLASS )
         mexErrMsgTxt("Input must be noncomplex single floating points.");
      }


   /* Check if EnableGPU is specified */
   enable_gpu = 1;

   /* Extract pointers to input matrices */
   float *image_ptr  = (float *) (mxGetData(prhs[0]));
   float *nodes_ptr  = (float *) (mxGetData(prhs[1]));

   /* Allocate out matrix */   
   mwSize mw_outimage_size[3];
   mw_outimage_size[0] = (mwSize)image_size[0];
   mw_outimage_size[1] = (mwSize)image_size[1];
   mw_outimage_size[2] = (mwSize)image_size[2];

   plhs[0] =  mxCreateNumericArray(3, mw_outimage_size, mxSINGLE_CLASS, mxREAL);
   float *outimage_ptr = (float *)(mxGetData(plhs[0]));

   /* Compute NMI gradient with respect of position of nodes of the splines */
   status = reg_array_resample_spline(outimage_ptr, image_ptr, nodes_ptr, image_size, nodes_size, spacing, enable_gpu);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while computing NMI gradient with respect of control point positions.");
   return;
}


