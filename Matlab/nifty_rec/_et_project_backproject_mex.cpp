
#include <mex.h>
#include "_et.h"

/*###############################################################################################*/
/* Matlab extensions */
/*###############################################################################################*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* Check for proper number of arguments. */
   /* Are expected reference, floating and bin number */ 
   if (nrhs!=3 ){
      mexErrMsgTxt("3 inputs required: Activity, Cameras, Point Spread Function");
   }

   const mxClassID classID = mxGetClassID(prhs[0]);
   const int dim=mxGetNumberOfDimensions(prhs[0]); 

   //...
   
   return;
}


