
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nifti1_io.h"


extern "C" void convolutionRow(float *h_Dst,float *h_Src,float *h_Kernel,int imageW,int imageH,int kernelR);
extern "C" void convolutionColumn(float *h_Dst,float *h_Src,float *h_Kernel,int imageW,int imageH,int kernelR);

int et_convolveSeparable2D(nifti_image* inputImage, float *kernelSeparated, int kernelLengthH, int kernelLengthW, nifti_image *outputImage, float background);
