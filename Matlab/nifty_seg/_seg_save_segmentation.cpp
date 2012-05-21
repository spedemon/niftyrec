#include <mex.h>
#include "_seg_array_interface.h"

#include <limits>
#include <string.h>
#include <math.h>
#include <cmath>

/*###############################################################################################*/
/* Matlab extensions */
/*###############################################################################################*/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   int status;
   char *filename;

   /* Check for proper number of arguments. */
   if (!(nrhs==1)){
      mexErrMsgTxt("1 input required: name for the segmentation nifti file");
   }
   /* Call NiftySeg function */
   filename = (char*) mxArrayToString(prhs[0]);
   status = seg_array_save_segmentation(filename);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while saving segmentation to file.");
   return;
}

