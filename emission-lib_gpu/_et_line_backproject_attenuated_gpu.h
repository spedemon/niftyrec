#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cutil.h>
#include "nifti1_io.h"

#define BLOCK 256

void et_line_backproject_attenuated_gpu(float **d_sinogram, float **d_backprojection, float **d_attenuation, int cam, nifti_image *backprojection);
