#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cutil.h>
#include "nifti1_io.h"

#define BLOCK 256

/*  d_B = d_B + d_A  */
void et_accumulate_gpu(float **d_A, float **d_B, nifti_image *n_image);
