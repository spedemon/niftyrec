#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cutil.h>
#include "nifti1_io.h"

void et_clear_accumulator_gpu(float **d_accumulator, nifti_image *accumulator);
