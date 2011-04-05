
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cutil.h>
#include "nifti1_io.h"

void tt_perspective_positionField_gpu(nifti_image *attenuation, float *image_origin, float *detector_origin, float *detector_size, float4 **positionField);
