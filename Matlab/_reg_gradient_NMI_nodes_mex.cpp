
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
   if (nrhs!=5){
      mexErrMsgTxt("5 inputs required: Target Image, Source Image, Spline nodes, Initial distance between nodes, Number of histogram bins");
   }

   const mxClassID cid_target  = mxGetClassID(prhs[0]);
   const int       dim_target  = mxGetNumberOfDimensions(prhs[0]); 

   const mxClassID cid_source   = mxGetClassID(prhs[1]);
   const int       dim_source   = mxGetNumberOfDimensions(prhs[1]); 

   int image_size[3];     //
   float spacing[3];      //
   int nodes_size[3];
   int enable_gpu;        // Flag for GPU Acceleration: 1 to enable, 0 to disable.

   int status = 1;        // Return 0 if everython ok, 1 if errors.

   image_size[0] = mxGetDimensions(prhs[0])[0];
   image_size[1] = mxGetDimensions(prhs[0])[1];
   image_size[2] = mxGetDimensions(prhs[0])[2];

   nodes_size[0] = mxGetDimensions(prhs[2])[0];
   nodes_size[1] = mxGetDimensions(prhs[2])[1];
   nodes_size[2] = mxGetDimensions(prhs[2])[2];

   spacing[0] = ((float*) mxGetData(prhs[3]))[0];
   spacing[1] = ((float*) mxGetData(prhs[3]))[1];
   spacing[2] = ((float*) mxGetData(prhs[3]))[2];

   /* The inputs must be noncomplex floating point matrices */
   for (int n=0; n<2; n++)
      {
      if ( mxGetClassID(prhs[n]) != mxSINGLE_CLASS )
         mexErrMsgTxt("Input must be noncomplex single floating points.");
      }


   /* Check if EnableGPU is specified */
   enable_gpu = 1;

   /* Extract pointers to input matrices */
   float *target_ptr  = (float *) (mxGetData(prhs[0]));
   float *source_ptr  = (float *) (mxGetData(prhs[1]));
   float *control_points_ptr = (float *) (mxGetData(prhs[2]));
   int binning = floor(((double*)(mxGetData(prhs[4])))[0]);

//nodes_size[0] = 128; 
//nodes_size[1] = 128; 
//nodes_size[2] = 128; 

   /* Allocate out matrix */   
/*   mwSize mw_gradient_size[4];
   mw_gradient_size[0] = (mwSize)nodes_size[0];
   mw_gradient_size[1] = (mwSize)nodes_size[1];
   mw_gradient_size[2] = (mwSize)nodes_size[2];
   mw_gradient_size[3] = 3;

   plhs[0] =  mxCreateNumericArray(4, mw_gradient_size, mxSINGLE_CLASS, mxREAL);
   float *gradient_ptr = (float *)(mxGetData(plhs[0]));
*/

fprintf(stderr, "- Nodes size: %d %d %d\n", nodes_size[0], nodes_size[1], nodes_size[2]);
fprintf(stderr, "- Image size: %d %d %d\n", image_size[0], image_size[1], image_size[2]);

mwSize mw_gradient_size[4];
mw_gradient_size[0] = 4;
mw_gradient_size[1] = (mwSize)nodes_size[0];
mw_gradient_size[2] = (mwSize)nodes_size[1];
mw_gradient_size[3] = (mwSize)nodes_size[2];

//mw_gradient_size[0] = 4;
plhs[0] =  mxCreateNumericArray(4, mw_gradient_size, mxSINGLE_CLASS, mxREAL);
float *gradient_ptr = (float *)(mxGetData(plhs[0]));

   /* Compute NMI gradient with respect of position of nodes of the splines */
   status = reg_array_gradient_NMI_nodes(target_ptr ,source_ptr, gradient_ptr, image_size, control_points_ptr, spacing, binning, enable_gpu);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while computing NMI gradient with respect of control point positions.");
   return;
}


