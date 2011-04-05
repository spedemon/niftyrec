


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
   int use_gpu = 0;

   /* Check for correct number of arguments. */
   if (nrhs!=3){
      mexErrMsgTxt("3 inputs required: Gradient wrt voxel intensity [N,M,K], Control points [Y,W,Z,3], Control Points Spacing [y,w,z]");
   }

   //extract pointers to input matrices
   float* gradient_ptr        = (float *) (mxGetData(prhs[0]));
   float* control_points_ptr  = (float *) (mxGetData(prhs[1]));
   float* grid_spacing_ptr    = (float *) (mxGetData(prhs[2]));

   //calculate size of control points grid, in order to allocate output
   int image_size[3];
   image_size[0] = mxGetDimensions(prhs[0])[0];
   image_size[1] = mxGetDimensions(prhs[0])[1];
   image_size[2] = mxGetDimensions(prhs[0])[2];

   int cp_size[3];
   reg_array_compute_control_point_size(cp_size, image_size, grid_spacing_ptr);

   //allocate output matrix for control points gradient  
   mwSize mw_cp_size[4];
   mw_cp_size[0] = (mwSize)cp_size[0];
   mw_cp_size[1] = (mwSize)cp_size[1];
   mw_cp_size[2] = (mwSize)cp_size[2];
   mw_cp_size[3] = 3;

fprintf(stderr,"%d %d %d \n",cp_size[0],cp_size[1],cp_size[2]);

   plhs[0] =  mxCreateNumericArray(4, mw_cp_size, mxSINGLE_CLASS, mxREAL);
   float *cp_gradient_ptr = (float *)(mxGetData(plhs[0]));

   /* Use gpu */
   use_gpu = 1;

   /* Compute control points from affine */
   status = reg_array_gradient_voxel_to_nodes(cp_gradient_ptr, gradient_ptr, control_points_ptr, image_size, cp_size, grid_spacing_ptr, use_gpu);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while computing gradient wrt control points position.");
   return;
}


