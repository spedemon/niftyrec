
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nifti1_io.h"

void et_line_integral_attenuated(nifti_image *inputImage, nifti_image *attenuationImage, nifti_image *sinoImage, int cam);
