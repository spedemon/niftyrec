#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cutil.h>
#include "nifti1_io.h"

void et_line_integral_attenuated_gpu(float **d_activity, float **d_attenuation, float **d_sinogram, int cam, nifti_image *);
