
#include<stdio.h>
#include<stdarg.h>

#include "_reg_tools.h"
#include "_reg_resampling.h"
#include "_reg_affineTransformation.h"

//#define _VERBOSE

extern "C" int fprintf_verbose (const char *__restrict __format, ...);
int et_create_rotation_matrix(mat44 *transformationMatrix, float theta_x, float theta_y, float theta_z, float center_x, float center_y, float center_z);
