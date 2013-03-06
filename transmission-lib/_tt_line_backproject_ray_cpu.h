/*
 *  _tt_line_backproject_ray_cpu.h
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */

#include <_tt_common.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//FIXME complete the CPU version (no CUDA)
//#include <cutil_inline.h>
#include <cutil_math.h>
#include "_reg_blocksize_gpu.h"

#include <vector_types.h>

extern "C" int tt_line_backproject_ray_cpu(float *out_backprojection, float *current_projection, float *invViewMatrix, uint2 detectorPixels, float3 sourcePosition, uint3 volumeVoxels, float3 volumeSize, float t_step, int interpolation); 


