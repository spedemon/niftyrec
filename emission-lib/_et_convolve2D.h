
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nifti1_io.h"

int et_convolve2D(nifti_image *inputImage, nifti_image *kernelImage, nifti_image *outImage, float background); 
